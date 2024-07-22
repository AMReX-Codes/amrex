
#ifndef BL_INDEXTYPE_H
#define BL_INDEXTYPE_H
#include <AMReX_Config.H>

#include <AMReX_IntVect.H>
#include <AMReX_SPACE.H>
#include <AMReX_Tuple.H>

#include <iosfwd>

namespace amrex {

/**
* \brief Type for defining CellIndex so that all IndexTypeND with different dimensions
* have the same CellIndex type.
*/
struct CellIndexEnum {
    //! The cell index type: one of CELL or NODE.
    enum CellIndex { CELL = 0, NODE = 1 };
};

/**
* \brief Cell-Based or Node-Based Indices
*
* The class IndexTypeND defines an index as being cell based or node (edge)
* based in each of the dim directions. This class defines an
* enumerated type CellIndex to be either CELL or NODE; i.e. each of the
* dim dimensions must be either CELL or NODE.
*/
template<int dim>
class IndexTypeND : public CellIndexEnum
{
public:
    static_assert(1 <= dim && dim <= 31, "The number of dimensions of IndexTypeND must be positive"
        " and less than 32");

    //! The default constructor
    AMREX_GPU_HOST_DEVICE
    constexpr IndexTypeND () noexcept = default;
    //! Construct an IndexTypeND identical to an IntVectND.
    AMREX_GPU_HOST_DEVICE
    explicit IndexTypeND (const IntVectND<dim>& iv) noexcept {
        for (int i=0; i<dim; ++i) {
            itype |= (iv[i] ? 1U : 0U) << i;
        }
    }
    /**
    * \brief Construct an IndexTypeND given an explicit CellIndex for
    * each direction.  The inputs for this constructor are N CellIndex,
    * where N is equal to the number of dimensions of the IndexTypeND.
    */
    template <class...Args,
        std::enable_if_t<
            (sizeof...(Args)+1 == dim) &&
            IsConvertible_v<CellIndex, Args...>,
        int> = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr IndexTypeND (CellIndex i, Args...js) noexcept {
        CellIndex locarr[dim] = {i, static_cast<CellIndex>(js)...};
        for (int s=0; s<dim; ++s) {
            itype |= ((locarr[s] == CellIndex::NODE) ? 1U : 0U) << s;
        }
    }
    // dtor, copy-ctor, copy-op=, move-ctor, and move-op= are compiler generated.

    //! Set IndexTypeND to be NODE based in direction dir.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set (int dir) noexcept { itype |= mask(dir); }
    //! Set IndexTypeND to be CELL based in direction dir.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void unset (int dir) noexcept { itype &= ~mask(dir); }
    //! True if IndexTypeND is NODE based in direction dir.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool test (int dir) const noexcept { return (itype & mask(dir)) != 0; }
    //! Set NODE based in all directions.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    void setall () noexcept { itype = (1U << dim) - 1; }
    //! Set CELL based in all directions.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void clear () noexcept { itype = 0; }
    //! True if this IndexTypeND is NODE based in any direction.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool any () const noexcept { return itype != 0; }
    //! True if IndexTypeND is valid.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool ok () const noexcept { return itype < (1U << dim); }
    //! Change from CELL to NODE or NODE to CELL in direction dir.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void flip (int i) noexcept { itype ^= mask(i); }
    //! True if IndexTypeNDs are identical.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator== (const IndexTypeND& t) const noexcept { return t.itype == itype; }
    //! True if IndexTypeNDs are not identical.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator!= (const IndexTypeND& t) const noexcept { return t.itype != itype; }
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator< (const IndexTypeND& t) const noexcept { return itype < t.itype; }
    //! True if the IndexTypeND is CELL based in all directions.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool cellCentered () const noexcept { return itype == 0; }
    //! True if the IndexTypeND is CELL based in dir-direction.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool cellCentered (int dir) const noexcept { return (itype & mask(dir)) == 0; }
    //! True if the IndexTypeND is NODE based in all directions.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool nodeCentered () const noexcept { return itype == (1U<<dim)-1; }
    //! True if the IndexTypeND is NODE based in dir-direction.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool nodeCentered (int dir) const noexcept { return (itype & mask(dir)) != 0; }
    //! Set IndexTypeND to CellIndex type t in direction dir.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setType (int dir, CellIndex t) noexcept { t == CELL ? unset(dir) : set(dir); }
    //! Returns  the CellIndex in direction dir.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    CellIndex ixType (int dir) const noexcept { return (CellIndex) ((itype & (1U<<dir)) >> dir); }
    //! Return an integer representing the IndexTypeND in direction dir.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator[] (int dir) const noexcept { return test(dir); }
    //! Returns the i'th CellIndex of the IndexTypeND. Used by structured bindings.
    template<std::size_t i>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    CellIndex get () const noexcept { static_assert(0<=i && i<dim); return ixType(i); }
    //! Fill an IntVectND of size dim with IndexTypeNDs.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> ixType () const noexcept {
        IntVectND<dim> retval(0);
        for (int i=0; i<dim; ++i) {
            retval[i] = test(i);
        }
        return retval;
    }
    //! Fill an IntVectND of size dim with IndexTypeNDs.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> toIntVect () const noexcept {
        IntVectND<dim> retval(0);
        for (int i=0; i<dim; ++i) {
            retval[i] = test(i);
        }
        return retval;
    }
    /**
    * \brief This static member function returns an IndexTypeND object of value
    * IndexTypeND::CELL.  It is provided as a convenience to our users
    * when defining a Box all of whose faces should be of type
    * IndexTypeND::CELL.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr IndexTypeND<dim> TheCellType () noexcept {
        return IndexTypeND<dim>{};
    }
    /**
    * \brief This static member function returns an IndexTypeND object of value
    * IndexTypeND::NODE.  It is provided as a convenience to our users
    * when defining a Box all of whose faces should be of type
    * IndexTypeND::NODE.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr IndexTypeND<dim> TheNodeType () noexcept {
        IndexTypeND<dim> retval{};
        retval.setall();
        return retval;
    }

    //! Return the size of this IndexTypeND.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr std::size_t size () noexcept {
        return static_cast<std::size_t>(dim);
    }

    //! Return the size of this IndexTypeND.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr int isize () noexcept {
        return dim;
    }

    using value_type = CellIndex;

    /**
    * \brief Returns a new IndexTypeND of size new_dim and
    * assigns the first new_dim values of this IndexTypeND to it.
    */
    template<int new_dim>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IndexTypeND<new_dim> shrink () const noexcept {
        static_assert(new_dim <= dim);
        IndexTypeND<new_dim> retval{};
        retval.getBits() = itype & ((1U << new_dim) - 1);
        return retval;
    }

    /**
    * \brief Returns a new IndexTypeND of size new_dim and
    * assigns all values of this IndexTypeND to it and fill_extra to the remaining elements.
    */
    template<int new_dim>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IndexTypeND<new_dim> expand (CellIndex fill_extra=CellIndex::CELL) const noexcept {
        static_assert(new_dim >= dim);
        IndexTypeND<new_dim> retval{};
        retval.getBits() = itype;
        if (fill_extra == CellIndex::NODE) {
            retval.getBits() |= (1U << new_dim) - (1U << dim);
        }
        return retval;
    }

    /**
    * \brief Returns a new IndexTypeND of size new_dim
    * by either shrinking or expanding this IndexTypeND
    */
    template<int new_dim>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IndexTypeND<new_dim> resize (CellIndex fill_extra=CellIndex::CELL) const noexcept {
        if constexpr (new_dim > dim) {
            return expand<new_dim>(fill_extra);
        } else {
            return shrink<new_dim>();
        }
    }

    //! Return the bit field representing the underlying data
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    unsigned int& getBits () noexcept { return itype; }

    //! Return the bit field representing the underlying data
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    const unsigned int& getBits () const noexcept { return itype; }

private:
    //! Returns 1<<k.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr unsigned int mask (int k) noexcept { return 1U<<k; }
    //! An integer holding the CellIndex in bits 0 - dim-1.
    unsigned int itype{0};
};

using IndexType = IndexTypeND<AMREX_SPACEDIM>;

// Template deduction guide for IndexTypeND
template<int dim>
AMREX_GPU_HOST_DEVICE // __device__ for HIP
IndexTypeND(const IntVectND<dim>&) -> IndexTypeND<dim>;

// Template deduction guide for IndexTypeND
template <class...Args,
    std::enable_if_t<
        IsConvertible_v<IndexType::CellIndex, Args...>,
    int> = 0>
AMREX_GPU_HOST_DEVICE // __device__ for HIP
IndexTypeND(IndexType::CellIndex, Args...) -> IndexTypeND<sizeof...(Args)+1>;

namespace detail {
    std::ostream& index_type_write (std::ostream& os, const unsigned int& iv, int dim);
    std::istream& index_type_read (std::istream& is, unsigned int& iv, int dim);

    template<class T, std::size_t...Ns>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T IndexTypeSplit_imp (T& retval, std::index_sequence<Ns...>, unsigned int src) noexcept {
        int dim_shift = 0;
        (
            (
                amrex::get<Ns>(retval).getBits() =
                    (src >> dim_shift) & ((1U << amrex::get<Ns>(retval).isize()) - 1),
                dim_shift += amrex::get<Ns>(retval).isize()
            ), ...
        );
        return retval;
    }
}

//! Write an IndexTypeND to an ostream in ASCII.
template<int dim>
std::ostream& operator<< (std::ostream& os, const IndexTypeND<dim>& it) {
    return detail::index_type_write(os, it.getBits(), dim);
}
//! Read an IndexTypeND from an istream.
template<int dim>
std::istream& operator>> (std::istream& is, IndexTypeND<dim>& it) {
    return detail::index_type_read(is, it.getBits(), dim);
}

/**
* \brief Returns a IndexTypeND obtained by concatenating the input IndexTypeNDs.
* The dimension of the return value equals the sum of the dimensions of the inputted IndexTypeNDs.
*/
template<int d, int...dims>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr IndexTypeND<detail::get_sum<d, dims...>()>
IndexTypeCat (const IndexTypeND<d>& v, const IndexTypeND<dims>&...vects) noexcept {
    IndexTypeND<detail::get_sum<d, dims...>()> retval{};
    retval.getBits() |= v.getBits();
    int dim_shift = v.isize();
    (
        (
            retval.getBits() |= (vects.getBits() << dim_shift),
            dim_shift += vects.isize()
        ), ...
    );
    return retval;
}

/**
* \brief Returns a tuple of IndexTypeND obtained by splitting the input IndexTypeND
* according to the dimensions specified by the template arguments.
*/
template<int d, int...dims>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr GpuTuple<IndexTypeND<d>, IndexTypeND<dims>...>
IndexTypeSplit (const IndexTypeND<detail::get_sum<d, dims...>()>& v) noexcept {
    GpuTuple<IndexTypeND<d>, IndexTypeND<dims>...> retval{};
    return detail::IndexTypeSplit_imp(retval,
                                      std::make_index_sequence<1 + sizeof...(dims)>(),
                                      v.getBits());
}

/**
* \brief Returns a new IndexTypeND of size new_dim and
* assigns the first new_dim values of v to it.
*/
template<int new_dim, int old_dim>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr IndexTypeND<new_dim>
IndexTypeShrink (const IndexTypeND<old_dim>& v) noexcept {
    return v.template shrink<new_dim>();
}

/**
* \brief Returns a new IndexTypeND of size new_dim and
* assigns all values of iv to it and fill_extra to the remaining elements.
*/
template<int new_dim, int old_dim>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr IndexTypeND<new_dim>
IndexTypeExpand (const IndexTypeND<old_dim>& v,
                 IndexType::CellIndex fill_extra=IndexType::CellIndex::CELL) noexcept {
    return v.template expand<new_dim>(fill_extra);
}

/**
* \brief Returns a new IndexTypeND of size new_dim
* by either shrinking or expanding iv.
*/
template<int new_dim, int old_dim>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr IndexTypeND<new_dim>
IndexTypeResize (const IndexTypeND<old_dim>& v,
                 IndexType::CellIndex fill_extra=IndexType::CellIndex::CELL) noexcept {
    return v.template resize<new_dim>(fill_extra);
}

} // namespace amrex

// Spcialize std::tuple_size for IndexTypeND. Used by structured bindings.
template<int dim>
struct std::tuple_size<amrex::IndexTypeND<dim>> {
    static constexpr std::size_t value = dim;
};

// Spcialize std::tuple_element for IndexTypeND. Used by structured bindings.
template<std::size_t s, int dim>
struct std::tuple_element<s, amrex::IndexTypeND<dim>> {
    using type = typename amrex::IndexTypeND<dim>::CellIndex;
};

#endif /*BL_INDEXTYPE_H*/
