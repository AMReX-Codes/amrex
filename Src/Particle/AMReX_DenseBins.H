#ifndef AMREX_DENSEBINS_H_
#define AMREX_DENSEBINS_H_

#include <AMReX_Gpu.H>
#include <AMReX_IntVect.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BinIterator.H>

namespace amrex
{

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

/**
 * \brief A container for storing items in a set of bins.
 *
 * The underlying data structure is an array of size nitems defining a
 * permutation of the items in the container that puts them in bin-sorted order,
 * plus an array of size nbins+1 that stores the offsets into the permutation
 * array where each bin starts.
 *
 * The storage for the bins is "dense" in the sense that users pass in
 * a Box that defines the space over which the bins will be defined, and
 * empty bins will still take up space.
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
        BL_PROFILE("DenseBins<T>::build");

        m_items = v;

        m_cells.resize(nitems);
        m_perm.resize(nitems);

        auto nbins = bx.numPts();
        m_counts.resize(0);
        m_counts.resize(nbins+1, 0);

        m_offsets.resize(0);
        m_offsets.resize(nbins+1);

        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        index_type* pcell   = m_cells.dataPtr();
        index_type* pcount  = m_counts.dataPtr();
        amrex::ParallelFor(nitems, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            auto iv = f(v[i]);
            auto iv3 = iv.dim3();
            int nx = hi.x-lo.x+1;
            int ny = hi.y-lo.y+1;
            int nz = hi.z-lo.z+1;
            index_type uix = amrex::min(nx-1,amrex::max(0,iv3.x));
            index_type uiy = amrex::min(ny-1,amrex::max(0,iv3.y));
            index_type uiz = amrex::min(nz-1,amrex::max(0,iv3.z));
            pcell[i] = (uix * ny + uiy) * nz + uiz;
            Gpu::Atomic::Add(&pcount[pcell[i]], index_type{ 1 });
        });

        Gpu::exclusive_scan(m_counts.begin(), m_counts.end(), m_offsets.begin());

        Gpu::copy(Gpu::deviceToDevice, m_offsets.begin(), m_offsets.end(), m_counts.begin());

        index_type* pperm = m_perm.dataPtr();
        constexpr index_type max_index = std::numeric_limits<index_type>::max();
        amrex::ParallelFor(nitems, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            index_type index = Gpu::Atomic::Inc(&pcount[pcell[i]], max_index);
            pperm[index] = i;
        });

        Gpu::Device::streamSynchronize();
    }

    template <typename N, typename F>
    void build (N nitems, T const* v, int nbins, F&& f)
    {
        BL_PROFILE("DenseBins<T>::build");

        m_items = v;

        m_cells.resize(nitems);
        m_perm.resize(nitems);

        m_counts.resize(0);
        m_counts.resize(nbins+1, 0);

        m_offsets.resize(0);
        m_offsets.resize(nbins+1);

        index_type* pcell   = m_cells.dataPtr();
        index_type* pcount  = m_counts.dataPtr();
        amrex::ParallelFor(nitems, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            pcell[i] = f(v[i]);
            Gpu::Atomic::Add(&pcount[pcell[i]], index_type{ 1 });
        });

        Gpu::exclusive_scan(m_counts.begin(), m_counts.end(), m_offsets.begin());

        Gpu::copy(Gpu::deviceToDevice, m_offsets.begin(), m_offsets.end(), m_counts.begin());

        index_type* pperm = m_perm.dataPtr();
        constexpr index_type max_index = std::numeric_limits<index_type>::max();
        amrex::ParallelFor(nitems, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            index_type index = Gpu::Atomic::Inc(&pcount[pcell[i]], max_index);
            pperm[index] = i;
        });

        Gpu::Device::streamSynchronize();
    }

    //! \brief the number of items in the container
    Long numItems () const noexcept { return m_perm.size(); }

    //! \brief the number of bins in the container
    Long numBins () const noexcept { return m_offsets.size()-1; }

    //! \brief returns the pointer to the permutation array
    index_type* permutationPtr () noexcept { return m_perm.dataPtr(); }

    //! \brief returns the pointer to the offsets array
    index_type* offsetsPtr () noexcept { return m_offsets.dataPtr(); }

    //! \brief returns const pointer to the permutation array
    const index_type* permutationPtr () const noexcept { return m_perm.dataPtr(); }

    //! \brief returns const pointer to the offsets array
    const index_type* offsetsPtr () const noexcept { return m_offsets.dataPtr(); }

    //! \brief returns a GPU-capable object that can create iterators over the items in a bin.
    DenseBinIteratorFactory<T> getBinIteratorFactory() const noexcept
    {
        return DenseBinIteratorFactory<T>(m_offsets, m_perm, m_items);
    }

private:

    const T* m_items;

    Gpu::DeviceVector<index_type> m_cells;
    Gpu::DeviceVector<index_type> m_counts;
    Gpu::DeviceVector<index_type> m_offsets;
    Gpu::DeviceVector<index_type> m_perm;
};

}

#endif
