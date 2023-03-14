#ifndef AMREX_DENSEBINS_H_
#define AMREX_DENSEBINS_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BinIterator.H>

#include <limits>

namespace amrex
{
namespace BinPolicy
{
    struct GPUBinPolicy {};
    struct OpenMPBinPolicy {};
    struct SerialBinPolicy {};

    static constexpr GPUBinPolicy      GPU{};
    static constexpr OpenMPBinPolicy   OpenMP{};
    static constexpr SerialBinPolicy   Serial{};

#ifdef AMREX_USE_GPU
    static constexpr GPUBinPolicy      Default{};
#else
    static constexpr OpenMPBinPolicy   Default{};
#endif
}

template <typename T>
struct DenseBinIteratorFactory
{
    using index_type = unsigned int;

    DenseBinIteratorFactory (const Gpu::DeviceVector<index_type>& offsets,
                             const Gpu::DeviceVector<index_type>& permutation,
                             const T* items)
        : m_offsets_ptr(offsets.dataPtr()),
          m_permutation_ptr(permutation.dataPtr()),
          m_items(items)
    {}

    AMREX_GPU_HOST_DEVICE
    BinIterator<T> getBinIterator(const int bin_number) const noexcept
    {
        return BinIterator<T>(bin_number, m_offsets_ptr, m_permutation_ptr, m_items);
    }

    const index_type* m_offsets_ptr;
    const index_type* m_permutation_ptr;
    const T* m_items;
};

namespace
{
    template <typename F, typename T, typename I>
    AMREX_GPU_HOST_DEVICE
    auto call_f(F const& f, T* pstruct_ptr, I& index)
        noexcept -> decltype(f(pstruct_ptr,index))
    {
        return f(pstruct_ptr[index],index);
    }

    template <typename F, typename T, typename I>
    AMREX_GPU_HOST_DEVICE
    auto call_f(F const& f, T* pstruct_ptr, I& index)
        noexcept -> decltype(f(pstruct_ptr[index]))
    {
        return f(pstruct_ptr[index]);
    }
}

/**
 * \brief A container for storing items in a set of bins.
 *
 * The underlying data structure is an array of size nitems defining a
 * permutation of the items in the container that puts them in bin-sorted order,
 * plus an array of size nbins+1 that stores the offsets into the permutation
 * array where each bin starts.
 *
 * The storage for the bins is "dense" in the sense that users pass in
 * either a Box or an integer that defines the space over which the bins will be defined,
 * and empty bins will still take up space.
 *
 * \tparam The type of items we hold
 *
 */
template <typename T>
class DenseBins
{
public:

    using BinIteratorFactory = DenseBinIteratorFactory<T>;
    using index_type = unsigned int;

    /**
     * \brief Populate the bins with a set of items.
     *
     * The algorithm is similar to a counting sort. First, we count the number
     * of items in each bin. Then, we perform a prefix sum on the resulting counts.
     * Finally, the set of partial sums is incremented in parallel using atomicInc,
     * which results in a permutation array that places the items in bin-sorted order.
     *
     * This version uses a 3D box to specificy the index space over which the bins
     * are defined.
     *
     * This overload uses the "default" parallelization strategy. If AMReX has been
     * compiled for GPU, it runs on the GPU. Otherwise, it executes in serial.
     *
     * \tparam N the 'size' type that can enumerate all the items
     * \tparam F a function that maps items to IntVect bins
     *
     * \param nitems the number of items to put in the bins
     * \param v pointer to the start of the items
     * \param bx the Box that defines the space over which the bins will be defined
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (N nitems, T const* v, const Box& bx, F&& f)
    {
        build(BinPolicy::Default, nitems, v, bx, std::forward<F>(f));
    }

    /**
     * \brief Populate the bins with a set of items.
     *
     * The algorithm is similar to a counting sort. First, we count the number
     * of items in each bin. Then, we perform a prefix sum on the resulting counts.
     * Finally, the set of partial sums is incremented in parallel using atomicInc,
     * which results in a permutation array that places the items in bin-sorted order.
     *
     * This version uses a 1D index space for the set of bins.
     *
     * This overload uses the "default" parallelization strategy. If AMReX has been
     * compiled for GPU, it runs on the GPU. Otherwise, it executes in serial.
     *
     * \tparam N the 'size' type that can enumerate all the items
     * \tparam F a function that maps items to IntVect bins
     *
     * \param nitems the number of items to put in the bins
     * \param v pointer to the start of the items
     * \param nbins the number of bins to use
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (N nitems, T const* v, int nbins, F&& f)
    {
        build(BinPolicy::Default, nitems, v, nbins, std::forward<F>(f));
    }

    /**
     * \brief Populate the bins with a set of items.
     *
     * The algorithm is similar to a counting sort. First, we count the number
     * of items in each bin. Then, we perform a prefix sum on the resulting counts.
     * Finally, the set of partial sums is incremented in parallel using atomicInc,
     * which results in a permutation array that places the items in bin-sorted order.
     *
     * This version uses a 3D box to specificy the index space over which the bins
     * are defined.
     *
     * This overload uses the "GPU" parallelization strategy. If AMReX has been
     * compiled for GPU, it runs on the GPU. Otherwise, it executes in serial.
     *
     * \tparam N the 'size' type that can enumerate all the items
     * \tparam F a function that maps items to IntVect bins
     *
     * \param nitems the number of items to put in the bins
     * \param v pointer to the start of the items
     * \param bx the Box that defines the space over which the bins will be defined
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (BinPolicy::GPUBinPolicy, N nitems, T const* v, const Box& bx, F&& f)
    {
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        build(BinPolicy::GPU, nitems, v, bx.numPts(),
              [=] AMREX_GPU_DEVICE (const T& t) noexcept
              {
                  auto iv = f(t);
                  auto iv3 = iv.dim3();
                  int nx = hi.x-lo.x+1;
                  int ny = hi.y-lo.y+1;
                  int nz = hi.z-lo.z+1;
                  index_type uix = amrex::min(nx-1,amrex::max(0,iv3.x));
                  index_type uiy = amrex::min(ny-1,amrex::max(0,iv3.y));
                  index_type uiz = amrex::min(nz-1,amrex::max(0,iv3.z));
                  return (uix * ny + uiy) * nz + uiz;
              });
    }

    /**
     * \brief Populate the bins with a set of items.
     *
     * The algorithm is similar to a counting sort. First, we count the number
     * of items in each bin. Then, we perform a prefix sum on the resulting counts.
     * Finally, the set of partial sums is incremented in parallel using atomicInc,
     * which results in a permutation array that places the items in bin-sorted order.
     *
     * This version uses a 1D index space for the set of bins.
     *
     * This overload uses the "GPU" parallelization strategy. If AMReX has been
     * compiled for GPU, it runs on the GPU. Otherwise, it executes in serial.
     *
     * \tparam N the 'size' type that can enumerate all the items
     * \tparam F a function that maps items to IntVect bins
     *
     * \param nitems the number of items to put in the bins
     * \param v pointer to the start of the items
     * \param nbins the number of bins to use
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (BinPolicy::GPUBinPolicy, N nitems, T const* v, int nbins, F&& f)
    {
        BL_PROFILE("DenseBins<T>::buildGPU");

        m_items = v;

        m_bins.resize(nitems);
        m_perm.resize(nitems);
        m_local_offsets.resize(nitems);

        m_counts.resize(0);
        m_counts.resize(nbins+1, 0);

        m_offsets.resize(0);
        m_offsets.resize(nbins+1);

        index_type* pbins   = m_bins.dataPtr();
        index_type* pcount  = m_counts.dataPtr();
        index_type* plocal_offsets  = m_local_offsets.dataPtr();
        amrex::ParallelFor(nitems, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            pbins[i] = call_f(f,v,i);
            index_type off = Gpu::Atomic::Add(&pcount[pbins[i]], index_type{ 1 });
            plocal_offsets[i] = off;
        });

        Gpu::exclusive_scan(m_counts.begin(), m_counts.end(), m_offsets.begin());

        index_type* pperm = m_perm.dataPtr();
        index_type* poffsets = m_offsets.dataPtr();
        amrex::ParallelFor(nitems, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            index_type index = poffsets[pbins[i]] + plocal_offsets[i];
            pperm[index] = i;
        });

        Gpu::Device::streamSynchronize();
    }



    template <typename F>
    void build2 (index_type nitems, index_type nbins, F&& f)
    {
        BL_PROFILE("DenseBins<T>::build2");

        constexpr index_type gpu_block_size = 1024;
        constexpr index_type gpu_block_size_m1 = gpu_block_size - 1;
        constexpr index_type llist_guard = std::numeric_limits<index_type>::max();

        m_items = nullptr;
        m_bins.resize(0);
        m_counts.resize(0);
        m_local_offsets.resize(0);
        m_offsets.resize(0);
        m_perm.resize(0);

        nbins = (nbins + gpu_block_size_m1) / gpu_block_size * gpu_block_size;

        Gpu::DeviceVector<index_type> llist_start(nbins, llist_guard);
        Gpu::DeviceVector<index_type> llist_next(nitems);
        m_perm.resize(nitems);
        Gpu::DeviceScalar<index_type> global_idx(0);

        index_type* pllist_start = llist_start.dataPtr();
        index_type* pllist_next = llist_next.dataPtr();
        index_type* pperm = m_perm.dataPtr();
        index_type* pglobal_idx = global_idx.dataPtr();

        amrex::ParallelFor(nitems, [=] AMREX_GPU_DEVICE (index_type i) noexcept
        {
            i = nitems - i - 1;
            pllist_next[i] = Gpu::Atomic::Exch(pllist_start + f(i), i);
        });

        amrex::launch<gpu_block_size>(nbins / gpu_block_size, Gpu::gpuStream(),
            [=] AMREX_GPU_DEVICE () {
                __shared__ index_type sdata[gpu_block_size];
                index_type current_idx = pllist_start[threadIdx.x + gpu_block_size * blockIdx.x];

                while (true) {
                    sdata[threadIdx.x] = index_type(current_idx != llist_guard);
                    index_type x = 0;

                    // simple block wide prefix sum
                    for (index_type i = 1; i<gpu_block_size; i*=2) {
                        __syncthreads();
                        if (threadIdx.x >= i) {
                            x = sdata[threadIdx.x - i];
                        }
                        __syncthreads();
                        if (threadIdx.x >= i) {
                            sdata[threadIdx.x] += x;
                        }
                    }
                    __syncthreads();
                    if (sdata[gpu_block_size_m1] == 0) {
                        break;
                    }
                    __syncthreads();
                    if (threadIdx.x == gpu_block_size_m1) {
                        x = sdata[gpu_block_size_m1];
                        sdata[gpu_block_size_m1] = Gpu::Atomic::Add(pglobal_idx, x);
                    }
                    __syncthreads();
                    if (threadIdx.x < gpu_block_size_m1) {
                        sdata[threadIdx.x] += sdata[gpu_block_size_m1];
                    }
                    __syncthreads();
                    if (threadIdx.x == gpu_block_size_m1) {
                        sdata[gpu_block_size_m1] += x;
                    }
                    __syncthreads();

                    if (current_idx != llist_guard) {
                        pperm[sdata[threadIdx.x] - 1] = current_idx;
                        current_idx = pllist_next[current_idx];
                    }
                }
            });

        Gpu::Device::streamSynchronize();
    }



    void build2 (index_type nitems, T const* v, Box bx, Geometry geom, const IntVect idx_type)
    {
        const IntVect refine_vect = max(idx_type, IntVect(1)).min(IntVect(2));
        const IntVect type_vect = idx_type - idx_type / 2 * 2;

        geom.refine(refine_vect);

        Box domain = geom.Domain();

        bx.convert(type_vect);
        domain.convert(type_vect);

        const RealVect dxi(geom.InvCellSize());
        const RealVect pos_offset = Real(0.5) * (RealVect(geom.ProbLo()) + RealVect(geom.ProbHi())
            - RealVect(geom.CellSize()) * RealVect(domain.smallEnd() + domain.bigEnd()));

        const int ref_product = AMREX_D_TERM(refine_vect[0], * refine_vect[1], * refine_vect[2]);
        const IntVect ref_offset(AMREX_D_DECL(1, refine_vect[0], refine_vect[0] * refine_vect[1]));

        build2(nitems, bx.numPts() * ref_product,
            [=] AMREX_GPU_DEVICE (index_type idx) noexcept
            {
                IntVect iv = ((v[idx].pos() - pos_offset) * dxi).round();

                IntVect iv_coarse = iv / refine_vect;
                IntVect iv_remainder = iv - iv_coarse * refine_vect;

                iv_coarse = iv_coarse.max(bx.smallEnd());
                iv_coarse = iv_coarse.min(bx.bigEnd());
                return bx.index(iv_coarse) + bx.numPts() * (iv_remainder * ref_offset).sum();
            });
    }




    /**
     * \brief Populate the bins with a set of items.
     *
     * The algorithm is similar to a counting sort. First, we count the number
     * of items in each bin. Then, we perform a prefix sum on the resulting counts.
     * Finally, the set of partial sums is incremented in parallel using atomicInc,
     * which results in a permutation array that places the items in bin-sorted order.
     *
     * This version uses a 3D box to specificy the index space over which the bins
     * are defined.
     *
     * This overload uses the "OpenMP" parallelization strategy and always runs on the
     * host. If AMReX has been compiled with OpenMP support, the executation will be
     * parallelized, otherwise it will be serial.
     *
     * \tparam N the 'size' type that can enumerate all the items
     * \tparam F a function that maps items to IntVect bins
     *
     * \param nitems the number of items to put in the bins
     * \param v pointer to the start of the items
     * \param bx the Box that defines the space over which the bins will be defined
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (BinPolicy::OpenMPBinPolicy, N nitems, T const* v, const Box& bx, F&& f)
    {
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        build(BinPolicy::OpenMP, nitems, v, bx.numPts(),
              [=] (const T& t) noexcept
              {
                  auto iv = f(t);
                  auto iv3 = iv.dim3();
                  int nx = hi.x-lo.x+1;
                  int ny = hi.y-lo.y+1;
                  int nz = hi.z-lo.z+1;
                  index_type uix = amrex::min(nx-1,amrex::max(0,iv3.x));
                  index_type uiy = amrex::min(ny-1,amrex::max(0,iv3.y));
                  index_type uiz = amrex::min(nz-1,amrex::max(0,iv3.z));
                  return (uix * ny + uiy) * nz + uiz;
              });
    }

    /**
     * \brief Populate the bins with a set of items.
     *
     * The algorithm is similar to a counting sort. First, we count the number
     * of items in each bin. Then, we perform a prefix sum on the resulting counts.
     * Finally, the set of partial sums is incremented in parallel using atomicInc,
     * which results in a permutation array that places the items in bin-sorted order.
     *
     * This version uses a 1D index space for the set of bins.
     *
     * This overload uses the "OpenMP" parallelization strategy and always runs on the
     * host. If AMReX has been compiled with OpenMP support, the executation will be
     * parallelized, otherwise it will be serial.
     *
     * \tparam N the 'size' type that can enumerate all the items
     * \tparam F a function that maps items to IntVect bins
     *
     * \param nitems the number of items to put in the bins
     * \param v pointer to the start of the items
     * \param nbins the number of bins to use
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (BinPolicy::OpenMPBinPolicy, N nitems, T const* v, int nbins, F&& f)
    {
        BL_PROFILE("DenseBins<T>::buildOpenMP");

        m_items = v;

        m_bins.resize(nitems);
        m_perm.resize(nitems);

        int nchunks = OpenMP::get_max_threads();
        int chunksize = nitems / nchunks;
        auto counts = (index_type*)(The_Arena()->alloc(nchunks*nbins*sizeof(index_type)));
        for (int i = 0; i < nbins*nchunks; ++i) { counts[i] = 0;}

        m_counts.resize(0);
        m_counts.resize(nbins+1, 0);

        m_offsets.resize(0);
        m_offsets.resize(nbins+1);

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
        for (int j = 0; j < nchunks; ++j) {
            int istart = j*chunksize;
            int istop = (j == nchunks-1) ? nitems : (j+1)*chunksize;
            for (int i = istart; i < istop; ++i) {
                m_bins[i] = call_f(f,v,i);
                ++counts[nbins*j+m_bins[i]];
            }
        }

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < nbins; ++i) {
            index_type total = 0;
            for (int j = 0; j < nchunks; ++j) {
                auto tmp = counts[nbins*j+i];
                counts[nbins*j+i] = total;
                total += tmp;
            }
            m_counts[i] = total;
        }

        // note - this part has to be serial
        m_offsets[0] = 0;
        for (int i = 0; i < nbins; ++i) {m_offsets[i+1] = m_offsets[i] + m_counts[i];}

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < nbins; ++i) {
            for (int j = 0; j < nchunks; ++j) {
                counts[nbins*j+i] += m_offsets[i];
            }
        }

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
        for (int j = 0; j < nchunks; ++j) {
            int istart = j*chunksize;
            int istop = (j == nchunks-1) ? nitems : (j+1)*chunksize;
            for (int i = istart; i < istop; ++i) {
                auto bid = m_bins[i];
                m_perm[counts[nbins*j+bid]++] = i;
            }
        }

        The_Arena()->free(counts);
    }

    /**
     * \brief Populate the bins with a set of items.
     *
     * The algorithm is similar to a counting sort. First, we count the number
     * of items in each bin. Then, we perform a prefix sum on the resulting counts.
     * Finally, the set of partial sums is incremented in parallel using atomicInc,
     * which results in a permutation array that places the items in bin-sorted order.
     *
     * This version uses a 3D box to specificy the index space over which the bins
     * are defined.
     *
     * This overload uses the "Serial" parallelization strategy. It always runs in
     * serial on the host, regardless of whether AMREX_USE_GPU or AMREX_USE_OMP
     * has been defined.
     *
     * \tparam N the 'size' type that can enumerate all the items
     * \tparam F a function that maps items to IntVect bins
     *
     * \param nitems the number of items to put in the bins
     * \param v pointer to the start of the items
     * \param bx the Box that defines the space over which the bins will be defined
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (BinPolicy::SerialBinPolicy, N nitems, T const* v, const Box& bx, F&& f)
    {
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        build(BinPolicy::Serial, nitems, v, bx.numPts(),
              [=] (const T& t) noexcept
              {
                  auto iv = f(t);
                  auto iv3 = iv.dim3();
                  int nx = hi.x-lo.x+1;
                  int ny = hi.y-lo.y+1;
                  int nz = hi.z-lo.z+1;
                  index_type uix = amrex::min(nx-1,amrex::max(0,iv3.x));
                  index_type uiy = amrex::min(ny-1,amrex::max(0,iv3.y));
                  index_type uiz = amrex::min(nz-1,amrex::max(0,iv3.z));
                  return (uix * ny + uiy) * nz + uiz;
              });
    }

    /**
     * \brief Populate the bins with a set of items.
     *
     * The algorithm is similar to a counting sort. First, we count the number
     * of items in each bin. Then, we perform a prefix sum on the resulting counts.
     * Finally, the set of partial sums is incremented in parallel using atomicInc,
     * which results in a permutation array that places the items in bin-sorted order.
     *
     * This version uses a 1D index space for the set of bins.
     *
     * This overload uses the "Serial" parallelization strategy. It always runs in
     * serial on the host, regardless of whether AMREX_USE_GPU or AMREX_USE_OMP
     * has been defined.
     *
     * \tparam N the 'size' type that can enumerate all the items
     * \tparam F a function that maps items to IntVect bins
     *
     * \param nitems the number of items to put in the bins
     * \param v pointer to the start of the items
     * \param nbins the number of bins to use
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (BinPolicy::SerialBinPolicy, N nitems, T const* v, int nbins, F&& f)
    {
        BL_PROFILE("DenseBins<T>::buildSerial");

        m_items = v;

        m_bins.resize(nitems);
        m_perm.resize(nitems);

        m_counts.resize(0);
        m_counts.resize(nbins+1, 0);

        m_offsets.resize(0);
        m_offsets.resize(nbins+1);

        for (int i = 0; i < nitems; ++i) {
            m_bins[i] = call_f(f,v,i);
            ++m_counts[m_bins[i]];
        }

        Gpu::exclusive_scan(m_counts.begin(), m_counts.end(), m_offsets.begin());

        Gpu::copy(Gpu::deviceToDevice, m_offsets.begin(), m_offsets.end(), m_counts.begin());

        for (int i = 0; i < nitems; ++i) {
            index_type index = m_counts[m_bins[i]]++;
            m_perm[index] = i;
        }
    }

public:

    //! \brief the number of items in the container
    Long numItems () const noexcept { return m_perm.size(); }

    //! \brief the number of bins in the container
    Long numBins () const noexcept { return m_offsets.size()-1; }

    //! \brief returns the pointer to the permutation array
    index_type* permutationPtr () noexcept { return m_perm.dataPtr(); }

    //! \brief returns the pointer to the offsets array
    index_type* offsetsPtr () noexcept { return m_offsets.dataPtr(); }

    //! \brief returns the pointer to the bins array
    index_type* binsPtr () noexcept { return m_bins.dataPtr(); }

    //! \brief returns const pointer to the permutation array
    const index_type* permutationPtr () const noexcept { return m_perm.dataPtr(); }

    //! \brief returns const pointer to the offsets array
    const index_type* offsetsPtr () const noexcept { return m_offsets.dataPtr(); }

    //! \brief returns the const pointer to the bins array
    const index_type* binsPtr () const noexcept { return m_bins.dataPtr(); }

    //! \brief returns a GPU-capable object that can create iterators over the items in a bin.
    DenseBinIteratorFactory<T> getBinIteratorFactory() const noexcept
    {
        return DenseBinIteratorFactory<T>(m_offsets, m_perm, m_items);
    }

private:

    const T* m_items;

    Gpu::DeviceVector<index_type> m_bins;
    Gpu::DeviceVector<index_type> m_counts;
    Gpu::DeviceVector<index_type> m_local_offsets;
    Gpu::DeviceVector<index_type> m_offsets;
    Gpu::DeviceVector<index_type> m_perm;
};

}

#endif
