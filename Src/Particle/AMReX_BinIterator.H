#ifndef AMREX_BINITERATOR_H_
#define AMREX_BINITERATOR_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>

#include <utility>

namespace amrex
{

template< class T >
constexpr decltype(T::is_particle_tile_data) IsParticleTileData () { // NOLINT(readability-const-return-type)
    return T::is_particle_tile_data;
}

template< class T, class...Args >
constexpr bool IsParticleTileData (Args...) {
    return false;
}

template <typename T>
struct BinIterator
{
    using index_type = unsigned int;

    using const_pointer_type = std::conditional_t<IsParticleTileData<T>(),
        T,
        const T*
    >;

    using const_pointer_input_type = std::conditional_t<IsParticleTileData<T>(),
        const T&,
        const T*
    >;

    struct iterator
    {
        AMREX_GPU_HOST_DEVICE
        iterator (index_type start, index_type stop,
                  const index_type* a_perm, const_pointer_input_type a_items)
            : m_items(a_items), m_perm(a_perm), m_index(start), m_stop(stop)
        {}

        AMREX_GPU_HOST_DEVICE
        void operator++ () { ++m_index;; }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE
        bool operator!= (iterator const& /*rhs*/) const { return m_index < m_stop; }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE
        auto operator* () const
        {
            return std::make_pair(m_perm[m_index], m_items[m_perm[m_index]]);
        }

    private:
        const_pointer_type m_items;
        const index_type* m_perm;
        index_type m_index;
        index_type m_stop;
    };

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    iterator begin () const
    {
        if (m_i == m_not_found) {
            return iterator(0, 0, m_permutation_ptr, m_items);
        }
        return iterator(m_offsets_ptr[m_i], m_offsets_ptr[m_i+1], m_permutation_ptr, m_items);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    iterator end () const {
        if (m_i == m_not_found) {
            return iterator(0, 0, m_permutation_ptr, m_items);
        }
        return iterator(m_offsets_ptr[m_i+1], m_offsets_ptr[m_i+1], m_permutation_ptr, m_items);
    }

    AMREX_GPU_HOST_DEVICE
    BinIterator (index_type i, const index_type *offsets_ptr,
                 const index_type *permutation_ptr, const_pointer_input_type items)
        : m_i(i), m_offsets_ptr(offsets_ptr), m_permutation_ptr(permutation_ptr), m_items(items)
    {}

private:
    index_type m_i;
    const index_type * m_offsets_ptr;
    const index_type * m_permutation_ptr;
    const_pointer_type m_items;

    static constexpr index_type m_not_found = std::numeric_limits<index_type>::max();
};

}

#endif
