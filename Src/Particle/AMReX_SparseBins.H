#ifndef AMREX_SPARSEBINS_H_
#define AMREX_SPARSEBINS_H_

#include <AMReX_Gpu.H>
#include <AMReX_IntVect.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BinIterator.H>

namespace amrex
{
    
template <typename T>
struct SparseBinIteratorFactory
{

    using index_type = unsigned int;
    
    SparseBinIteratorFactory (const Gpu::DeviceVector<index_type>& bins,
                              const Gpu::DeviceVector<index_type>& offsets, 
                              const Gpu::DeviceVector<index_type>& permutation,
                              const T* items)
        : m_bins_ptr(bins.dataPtr()),
          m_offsets_ptr(offsets.dataPtr()),
          m_permutation_ptr(permutation.dataPtr()),
          m_items(items), m_num_bins(bins.size())
    {}

    AMREX_GPU_HOST_DEVICE
    index_type getIndex (const index_type bin_number) const noexcept
    {
        if (m_num_bins == 1) return 0;

        index_type lo = 0;
        index_type hi = m_num_bins - 1;

        while (lo < hi) {
            if (m_bins_ptr[lo] == bin_number) return lo;
            if (m_bins_ptr[hi] == bin_number) return hi;
            
            index_type mid = (lo + hi) / 2;
            index_type mid_value = m_bins_ptr[mid];
            if (mid_value == bin_number) return mid;
            
            mid_value < bin_number ? lo = mid+1 : hi = mid;            
        }
       
        return m_not_found;
    }    
    
    AMREX_GPU_HOST_DEVICE
    BinIterator<T> getBinIterator (const index_type bin_number) const noexcept
    {
        auto index = getIndex(bin_number);
        return BinIterator<T>(index, m_offsets_ptr, m_permutation_ptr, m_items);
    }

    const index_type* m_bins_ptr;
    const index_type* m_offsets_ptr;
    const index_type* m_permutation_ptr;
    const T* m_items;
    index_type m_num_bins;
    
    static constexpr index_type m_not_found = std::numeric_limits<index_type>::max(); 
};
    
/**
 * \brief A container for storing items in a set of bins using "sparse" storage.
 * 
 * The underlying data structure consists of three arrays. The first is a sorted
 * array of bin numbers with a size equal to the number of non-zero bins. The 
 * second is an array of size nitems defining a permutation of the items in the container 
 * that puts them in bin-sorted order. Finally, there is an offsets array that stores,
 * for each non-zero bin, the offsets into the permutation array where each bin starts.
 *
 * The storage for the bins is "sparse" in the sense that users pass in
 * a Box that defines the space over which the bins will be defined, and
 * empty bins will still take up space.
 *
 * \tparam The type of items we hold
 *
 */    
template <typename T>
class SparseBins
{
public:

    using BinIteratorFactory = SparseBinIteratorFactory<T>;
    using bin_type = IntVect;
    using index_type = unsigned int;

    /**
     * \brief Populate the bins with a set of items.
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
        BL_PROFILE("SparseBins<T>::build");

        m_items = v;
        
        Gpu::HostVector<index_type> host_cells(nitems);
        std::map<index_type, index_type> bins_map;
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        for (int i = 0; i < nitems; ++i)
        {
            bin_type iv = f(v[i]);
            auto iv3 = iv.dim3();
            int nx = hi.x-lo.x+1;
            int ny = hi.y-lo.y+1;
            int nz = hi.z-lo.z+1;
            index_type uix = amrex::min(nx-1,amrex::max(0,iv3.x));
            index_type uiy = amrex::min(ny-1,amrex::max(0,iv3.y));
            index_type uiz = amrex::min(nz-1,amrex::max(0,iv3.z));
            host_cells[i] = (uix * ny + uiy) * nz + uiz;
            bins_map[host_cells[i]] += 1;
        }

        Gpu::HostVector<index_type> host_perm(nitems);
        std::iota(host_perm.begin(), host_perm.end(), 0);
        std::sort(host_perm.begin(), host_perm.end(),
                  [&](int i, int j) {return host_cells[i] < host_cells[j];});
        
        Gpu::HostVector<index_type> host_bins;
        Gpu::HostVector<index_type> host_offsets = {0};
        for (const auto& kv : bins_map)
        {
            host_bins.push_back(kv.first);
            host_offsets.push_back(host_offsets.back() + kv.second);
        }
        
        m_bins.resize(host_bins.size());
        Gpu::copy(Gpu::hostToDevice, host_bins.begin(), host_bins.end(), m_bins.begin());
        
        m_offsets.resize(host_offsets.size());
        Gpu::copy(Gpu::hostToDevice, host_offsets.begin(), host_offsets.end(), m_offsets.begin());

        m_perm.resize(host_perm.size());
        Gpu::copy(Gpu::hostToDevice, host_perm.begin(), host_perm.end(), m_perm.begin());
    }

    //! \brief the number of items in the container
    Long numItems () const noexcept { return m_perm.size(); }

    //! \brief the number of bins in the container
    Long numBins () const noexcept { return m_offsets.size()-1; }

    //! \brief returns the pointer to the permutation array 
    index_type* permutationPtr () noexcept { return m_perm.dataPtr(); }

    //! \brief returns the pointer to the offsets array
    index_type* offsetsPtr () noexcept { return m_offsets.dataPtr(); }

    //! \brief returns the pointer to the array of non-zero bins
    index_type* getNonZeroBinsPtr() noexcept { return m_bins.dataPtr(); }
    
    //! \brief returns const pointer to the permutation array 
    const index_type* permutationPtr () const noexcept { return m_perm.dataPtr(); }

    //! \brief returns const pointer to the offsets array
    const index_type* offsetsPtr () const noexcept { return m_offsets.dataPtr(); }

    //! \brief returns the pointer to the array of non-zero bins
    const index_type* getNonZeroBinsPtr() const noexcept { return m_bins.dataPtr(); }

    //! \brief returns a GPU-capable object that can create iterators over the items in a bin.
    SparseBinIteratorFactory<T> getBinIteratorFactory() const noexcept
    {
        return SparseBinIteratorFactory<T>(m_bins, m_offsets, m_perm, m_items);
    }
    
private:

    const T* m_items;
    
    Gpu::DeviceVector<index_type> m_bins;
    Gpu::DeviceVector<index_type> m_offsets;
    Gpu::DeviceVector<index_type> m_perm;
};

}

#endif
