#ifndef AMREX_SPARSEBINS_H_
#define AMREX_SPARSEBINS_H_

#include <AMReX_Gpu.H>

#include <utility>

namespace amrex
{

template <typename T>
struct BinIterator
{
    struct iterator
    {
        AMREX_GPU_HOST_DEVICE
        iterator (int start, int stop, const unsigned int* a_perm, const T* a_items) 
            : m_index(start), m_stop(stop), m_perm(a_perm), m_items(a_items)   
        {}
        
        AMREX_GPU_HOST_DEVICE
        void operator++ () { ++m_index;; }
        
        AMREX_GPU_HOST_DEVICE
        bool operator!= (iterator const& rhs) const { return m_index < m_stop; }
        
        AMREX_GPU_HOST_DEVICE
        std::pair<unsigned int, T> operator* () const
        {
            return std::make_pair(m_perm[m_index], m_items[m_perm[m_index]]);
        }
        
    private:
        const T* m_items;
        const unsigned int* m_perm;
        int m_index;
        int m_stop;
    };
    
    AMREX_GPU_HOST_DEVICE
    iterator begin () const {
        return iterator(m_offsets_ptr[m_i], m_offsets_ptr[m_i+1], m_permutation_ptr, m_items);
    }

    AMREX_GPU_HOST_DEVICE
    iterator end () const { 
        return iterator(m_offsets_ptr[m_i+1], m_offsets_ptr[m_i+1], m_permutation_ptr, m_items);
    }

    AMREX_GPU_HOST_DEVICE    
    BinIterator (int i, const unsigned int *offsets_ptr, const unsigned int *permutation_ptr, const T* items)
        : m_i(i), m_offsets_ptr(offsets_ptr), m_permutation_ptr(permutation_ptr), m_items(items)
    {}
    
private:
    const int m_i;
    const unsigned int * m_offsets_ptr;
    const unsigned int * m_permutation_ptr;
    const T * m_items;
};

}

#endif
