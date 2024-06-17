
#ifndef BL_INDEXTYPE_H
#define BL_INDEXTYPE_H
#include <AMReX_Config.H>

#include <AMReX_IntVect.H>
#include <AMReX_SPACE.H>

#include <iosfwd>

namespace amrex {

/**
* \brief Cell-Based or Node-Based Indices
*
* The class IndexTypeND defines an index as being cell based or node (edge)
* based in each of the dim directions. This class defines an
* enumerated type CellIndex to be either CELL or NODE; i.e. each of the
* dim dimensions must be either CELL or NODE.
*/
template<int dim>
class IndexTypeND
{
public:
    static_assert(1 <= dim && dim <= 31, "The number of dimensions of IndexTypeND must be positive"
        " and less than 32");

    //! The cell index type: one of CELL or NODE.
    enum CellIndex { CELL = 0, NODE = 1 };
    //! The default constructor
    AMREX_GPU_HOST_DEVICE
    constexpr IndexTypeND () noexcept = default;
    //! Construct an IndexTypeND identical to an IntVectND.
    AMREX_GPU_HOST_DEVICE
    explicit IndexTypeND (const IntVectND<dim>& iv) noexcept {
        itype = 0;
        for (int i=0; i<dim; ++i) {
            itype |= (iv[i] ? 1u : 0u) << dim;
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
        itype = 0;
        for (int s=0; s<dim; ++s) {
            itype |= (locarr[s] ? 1u : 0u) << dim;
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
    void setall () noexcept { itype = (1u << dim) - 1; }
    //! Set CELL based in all directions.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void clear () noexcept { itype = 0; }
    //! True if this IndexTypeND is NODE based in any direction.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool any () const noexcept { return itype != 0; }
    //! True if IndexTypeND is valid.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool ok () const noexcept { return itype < (1u << dim); }
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
    bool nodeCentered () const noexcept { return itype == (1u<<dim)-1; }
    //! True if the IndexTypeND is NODE based in dir-direction.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool nodeCentered (int dir) const noexcept { return (itype & mask(dir)) != 0; }
    //! Set IndexTypeND to CellIndex type t in direction dir.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setType (int dir, CellIndex t) noexcept { t == CELL ? unset(dir) : set(dir); }
    //! Returns  the CellIndex in direction dir.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    CellIndex ixType (int dir) const noexcept { return (CellIndex) ((itype & (1u<<dir)) >> dir); }
    //! Return an integer representing the IndexTypeND in direction dir.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator[] (int dir) const noexcept { return test(dir); }
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

    template<int d> friend std::ostream& operator<< (std::ostream& os, const IndexTypeND<d>& itype);
    template<int d> friend std::istream& operator>> (std::istream& is, IndexTypeND<d>& itype);

private:
    //! Returns 1<<k.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr unsigned int mask (int k) noexcept { return 1u<<k; }
    //! An integer holding the CellIndex in bits 0 - dim-1.
    unsigned int itype{0};
};

using IndexType = IndexTypeND<AMREX_SPACEDIM>;

namespace detail {
    std::ostream& index_type_write (std::ostream& os, const unsigned int& iv, int dim);
    std::istream& index_type_read (std::istream& is, unsigned int& iv, int dim);
}

//! Write an IndexTypeND to an ostream in ASCII.
template<int dim>
std::ostream& operator<< (std::ostream& os, const IndexTypeND<dim>& itype) {
    return detail::index_type_write(os, itype.itype, dim);
}
//! Read an IndexTypeND from an istream.
template<int dim>
std::istream& operator>> (std::istream& is, IndexTypeND<dim>& itype) {
    return detail::index_type_read(is, itype.itype, dim);
}

}

#endif /*BL_INDEXTYPE_H*/
