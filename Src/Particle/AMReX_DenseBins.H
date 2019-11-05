#ifndef AMREX_DENSEBINS_H_
#define AMREX_DENSEBINS_H_

#include <AMReX_Gpu.H>
#include <AMReX_IntVect.H>

namespace amrex
{

/**
 * \brief A container for storing items in a set of bins.
 * 
 * The underlying data structure is an array of size nitems defining a
 * permutation of the items in the container that puts them in bin-sorted order,
 * plus an array of size nbins that stores the offsets into the permutation
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
    using bin_type = IntVect;
    using index_type = unsigned int;
    
public:

    /**
     * \brief Populate the bins with a set of items.
     *
     * The algorithm is similar to a counting sort. First, we count the number
     * of items in each bin. Then, we perform a prefix sum on the resulting counts.
     * Finally, the set of partial sums is incremented in parallel using atomicInc,
     * which results in a permutationa array that places the items in bin-sorted order.
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
    void build (N nitems, T const* v, const Box& bx, F f,
                const bin_type& bin_size = bin_type{ AMREX_D_DECL(1, 1, 1) })
    {
        BL_PROFILE("DenseBins<T>::build");

        m_cells.resize(nitems);
        m_perm.resize(nitems);
        
        m_counts.resize(0);
        m_counts.resize(bx.numPts(), 0);

        m_offsets.resize(0);
        m_offsets.resize(bx.numPts()+1);

        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        IntVect num_bins = (bx.bigEnd() - bx.smallEnd() + bin_size) / bin_size;
        auto nbins = num_bins.dim3();
        auto bsize = bin_size.dim3();

        // clamp to 1 for AMREX_SPACEDIM < 3
        nbins.x = amrex::max(nbins.x, 1);
        nbins.y = amrex::max(nbins.y, 1);
        nbins.z = amrex::max(nbins.z, 1);

        bsize.x = amrex::max(nbins.x, 1);
        bsize.y = amrex::max(nbins.y, 1);
        bsize.z = amrex::max(nbins.z, 1);
            
        index_type* pcell   = m_cells.dataPtr();
        index_type* pcount  = m_counts.dataPtr();
        index_type* poffset = m_offsets.dataPtr();
        AMREX_FOR_1D ( nitems, i,
        {
            bin_type iv = f(v[i]);
            auto ind = iv.dim3();
            int ix = ind.x / bsize.x;
            int iy = ind.y / bsize.y;
            int iz = ind.z / bsize.z;            
            int nx = nbins.x;
            int ny = nbins.y;
            int nz = nbins.z;
            index_type uix = min(nx-1,max(0,ix));
            index_type uiy = min(ny-1,max(0,iy));
            index_type uiz = min(nz-1,max(0,iz));
            pcell[i] = (uix * ny + uiy) * nz + uiz;
            Gpu::Atomic::Add(&pcount[pcell[i]], index_type{ 1 });
        });

        Gpu::exclusive_scan(m_counts.begin(), m_counts.end()+1, m_offsets.begin());
        
        Gpu::thrust_copy(m_offsets.begin(), m_offsets.end()-1, m_counts.begin());

        index_type* pperm = m_perm.dataPtr();
        constexpr index_type max_index = std::numeric_limits<index_type>::max();
        AMREX_FOR_1D ( nitems, i,
        {
            index_type index = Gpu::Atomic::Inc(&pcount[pcell[i]], max_index);
            pperm[index] = i;
        });
    }

    //! \brief the number of items in the container
    long numItems () const noexcept { return m_perm.size(); }

    //! \brief the number of bins in the container
    long numBins () const noexcept { return m_offsets.size(); }
    
    //! \brief returns the pointer to the permutation array 
    index_type* permutationPtr() noexcept { return m_perm.dataPtr(); }

    //! \brief returns the pointer to the offsets array
    index_type* offsetsPtr() noexcept { return m_offsets.dataPtr(); }

    //! \brief returns const pointer to the permutation array 
    const index_type* permutationPtr() const noexcept { return m_perm.dataPtr(); }

    //! \brief returns const pointer to the offsets array
    const index_type* offsetsPtr() const noexcept { return m_offsets.dataPtr(); }

    
private:

    Gpu::DeviceVector<index_type> m_cells;
    Gpu::DeviceVector<index_type> m_counts;
    Gpu::DeviceVector<index_type> m_offsets;
    Gpu::DeviceVector<index_type> m_perm;
};

}

#endif
