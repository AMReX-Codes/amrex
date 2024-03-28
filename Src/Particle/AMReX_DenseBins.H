#ifndef AMREX_DENSEBINS_H_
#define AMREX_DENSEBINS_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>
#include <AMReX_Scan.H>
#include <AMReX_IntVect.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BinIterator.H>

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
    using index_type = int;

    using const_pointer_type = std::conditional_t<IsParticleTileData<T>(),
        T,
        const T*
    >;

    DenseBinIteratorFactory (const Gpu::DeviceVector<index_type>& offsets,
                             const Gpu::DeviceVector<index_type>& permutation,
                             const T* items)
        : m_offsets_ptr(offsets.dataPtr()),
          m_permutation_ptr(permutation.dataPtr()),
          m_items(items)
    {}

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    BinIterator<T> getBinIterator(const int bin_number) const noexcept
    {
        return BinIterator<T>(bin_number, m_offsets_ptr, m_permutation_ptr, m_items);
    }

    const index_type* m_offsets_ptr;
    const index_type* m_permutation_ptr;
    const_pointer_type m_items;
};


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
    using index_type = int;

    using const_pointer_type = std::conditional_t<IsParticleTileData<T>(),
        T,
        const T*
    >;

    using const_pointer_input_type = std::conditional_t<IsParticleTileData<T>(),
        const T&,
        const T*
    >;

private:

    template <typename F, typename I>
    AMREX_GPU_HOST_DEVICE
    static auto call_f (F const& f, const_pointer_input_type v, I& index) {
        if constexpr (IsCallable<F, decltype(v), I>::value) {
            return f(v, index);
        } else {
            return f(v[index]);
        }
    }

public:

    /**
     * \brief Populate the bins with a set of items.
     *
     * The algorithm is similar to a counting sort. First, we count the number
     * of items in each bin. Then, we perform a prefix sum on the resulting counts.
     * Finally, the set of partial sums is incremented in parallel using atomicInc,
     * which results in a permutation array that places the items in bin-sorted order.
     *
     * This version uses a 3D box to specify the index space over which the bins
     * are defined.
     *
     * This overload uses the "default" parallelization strategy. If AMReX has been
     * compiled for GPU, it runs on the GPU. Otherwise, it executes in serial.
     *
     * \tparam N the 'size' type that can enumerate all the items
     * \tparam F a function that maps items to IntVect bins
     *
     * \param nitems the number of items to put in the bins
     * \param v ParticleTileData or pointer to the start of the items
     * \param bx the Box that defines the space over which the bins will be defined
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (N nitems, const_pointer_input_type v, const Box& bx, F&& f)
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
     * \param v ParticleTileData or pointer to the start of the items
     * \param nbins the number of bins to use
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (N nitems, const_pointer_input_type v, int nbins, F&& f)
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
     * This version uses a 3D box to specify the index space over which the bins
     * are defined.
     *
     * This overload uses the "GPU" parallelization strategy. If AMReX has been
     * compiled for GPU, it runs on the GPU. Otherwise, it executes in serial.
     *
     * \tparam N the 'size' type that can enumerate all the items
     * \tparam F a function that maps items to IntVect bins
     *
     * \param nitems the number of items to put in the bins
     * \param v ParticleTileData or pointer to the start of the items
     * \param bx the Box that defines the space over which the bins will be defined
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (BinPolicy::GPUBinPolicy, N nitems, const_pointer_input_type v, const Box& bx, F const& f)
    {
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        build(BinPolicy::GPU, nitems, v, bx.numPts(),
              [=] AMREX_GPU_DEVICE (const_pointer_type t, index_type i) noexcept
              {
                  auto iv = call_f(f,t,i);
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
     * \param v ParticleTileData or pointer to the start of the items
     * \param nbins the number of bins to use
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (BinPolicy::GPUBinPolicy, N nitems, const_pointer_input_type v, int nbins, F const& f)
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

    /**
     * \brief Populate the bins with a set of items.
     *
     * The algorithm is similar to a counting sort. First, we count the number
     * of items in each bin. Then, we perform a prefix sum on the resulting counts.
     * Finally, the set of partial sums is incremented in parallel using atomicInc,
     * which results in a permutation array that places the items in bin-sorted order.
     *
     * This version uses a 3D box to specify the index space over which the bins
     * are defined.
     *
     * This overload uses the "OpenMP" parallelization strategy and always runs on the
     * host. If AMReX has been compiled with OpenMP support, the execution will be
     * parallelized, otherwise it will be serial.
     *
     * \tparam N the 'size' type that can enumerate all the items
     * \tparam F a function that maps items to IntVect bins
     *
     * \param nitems the number of items to put in the bins
     * \param v ParticleTileData or pointer to the start of the items
     * \param bx the Box that defines the space over which the bins will be defined
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (BinPolicy::OpenMPBinPolicy, N nitems, const_pointer_input_type v, const Box& bx, F const& f)
    {
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        build(BinPolicy::OpenMP, nitems, v, bx.numPts(),
              [=] (const_pointer_type t, index_type i) noexcept
              {
                  auto iv = call_f(f,t,i);
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
     * host. If AMReX has been compiled with OpenMP support, the execution will be
     * parallelized, otherwise it will be serial.
     *
     * \tparam N the 'size' type that can enumerate all the items
     * \tparam F a function that maps items to IntVect bins
     *
     * \param nitems the number of items to put in the bins
     * \param v ParticleTileData or pointer to the start of the items
     * \param nbins the number of bins to use
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (BinPolicy::OpenMPBinPolicy, N nitems, const_pointer_input_type v, int nbins, F const& f)
    {
        BL_PROFILE("DenseBins<T>::buildOpenMP");

        m_items = v;

        m_bins.resize(nitems);
        m_perm.resize(nitems);

        int nchunks = OpenMP::get_max_threads();
        int chunksize = nitems / nchunks;
        auto* counts = (index_type*)(The_Arena()->alloc(nchunks*nbins*sizeof(index_type)));
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
     * This version uses a 3D box to specify the index space over which the bins
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
     * \param v ParticleTileData or pointer to the start of the items
     * \param bx the Box that defines the space over which the bins will be defined
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (BinPolicy::SerialBinPolicy, N nitems, const_pointer_input_type v, const Box& bx, F const& f)
    {
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        build(BinPolicy::Serial, nitems, v, bx.numPts(),
              [=] (const_pointer_type t, index_type i) noexcept
              {
                  auto iv = call_f(f,t,i);
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
     * \param v ParticleTileData or pointer to the start of the items
     * \param nbins the number of bins to use
     * \param f a function object that maps items to bins
     */
    template <typename N, typename F>
    void build (BinPolicy::SerialBinPolicy, N nitems, const_pointer_input_type v, int nbins, F const& f)
    {
        BL_PROFILE("DenseBins<T>::buildSerial");

        m_items = v;

        m_bins.resize(nitems);
        m_perm.resize(nitems);

        m_counts.resize(0);
        m_counts.resize(nbins+1, 0);

        m_offsets.resize(0);
        m_offsets.resize(nbins+1);

        for (N i = 0; i < nitems; ++i) {
            m_bins[i] = call_f(f,v,i);
            ++m_counts[m_bins[i]];
        }

        Gpu::exclusive_scan(m_counts.begin(), m_counts.end(), m_offsets.begin());

        Gpu::copy(Gpu::deviceToDevice, m_offsets.begin(), m_offsets.end(), m_counts.begin());

        for (N i = 0; i < nitems; ++i) {
            index_type index = m_counts[m_bins[i]]++;
            m_perm[index] = i;
        }
    }

    //! \brief the number of items in the container
    [[nodiscard]] Long numItems () const noexcept { return m_perm.size(); }

    //! \brief the number of bins in the container
    [[nodiscard]] Long numBins () const noexcept { return m_offsets.size()-1; }

    //! \brief returns the pointer to the permutation array
    [[nodiscard]] index_type* permutationPtr () noexcept { return m_perm.dataPtr(); }

    //! \brief returns the pointer to the offsets array
    [[nodiscard]] index_type* offsetsPtr () noexcept { return m_offsets.dataPtr(); }

    //! \brief returns the pointer to the bins array
    [[nodiscard]] index_type* binsPtr () noexcept { return m_bins.dataPtr(); }

    //! \brief returns const pointer to the permutation array
    [[nodiscard]] const index_type* permutationPtr () const noexcept { return m_perm.dataPtr(); }

    //! \brief returns const pointer to the offsets array
    [[nodiscard]] const index_type* offsetsPtr () const noexcept { return m_offsets.dataPtr(); }

    //! \brief returns the const pointer to the bins array
    [[nodiscard]] const index_type* binsPtr () const noexcept { return m_bins.dataPtr(); }

    //! \brief returns a GPU-capable object that can create iterators over the items in a bin.
    [[nodiscard]] DenseBinIteratorFactory<T> getBinIteratorFactory() const noexcept
    {
        return DenseBinIteratorFactory<T>(m_offsets, m_perm, m_items);
    }

private:

    const_pointer_type m_items;

    Gpu::DeviceVector<index_type> m_bins;
    Gpu::DeviceVector<index_type> m_counts;
    Gpu::DeviceVector<index_type> m_local_offsets;
    Gpu::DeviceVector<index_type> m_offsets;
    Gpu::DeviceVector<index_type> m_perm;
};

}

#endif
