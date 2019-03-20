
#ifndef BL_INDEXTYPE_H
#define BL_INDEXTYPE_H

#include <iosfwd>

#include <AMReX_ccse-mpi.H>
#include <AMReX_IntVect.H>
#include <AMReX_SPACE.H>

namespace amrex {

/**
* \brief Cell-Based or Node-Based Indices
*
* The class IndexType defines an index as being cell based or node (edge)
* based in each of the AMREX_SPACEDIM directions.  This class defines an
* enumerated type CellIndex to be either CELL or NODE; i.e. each of the
* AMREX_SPACEDIM dimensions must be either CELL or NODE.
*/

class IndexType
{
    friend MPI_Datatype ParallelDescriptor::Mpi_typemap<IndexType>::type();

public:
    //! The cell index type: one of CELL or NODE.
    enum CellIndex { CELL = 0, NODE = 1 };
    //! The default constructor
    AMREX_GPU_HOST_DEVICE
    constexpr IndexType () noexcept : itype(0) {}
    //! Construct an IndexType identical to an IntVect.
    AMREX_GPU_HOST_DEVICE
    explicit IndexType (const IntVect& iv) noexcept
        : itype(AMREX_D_TERM((iv[0]?1:0), | ((iv[1]?1:0)<<1), | ((iv[2]?1:0)<<2)))
        {}
    /**
    * \brief Construct an IndexType given an explicit CellIndex for
    * each direction.  AMREX_D_DECL is a macro that sets the constructor
    * to take AMREX_SPACEDIM arguments.
    */
    AMREX_GPU_HOST_DEVICE
    constexpr IndexType (AMREX_D_DECL(CellIndex i, CellIndex j, CellIndex k)) noexcept
        : itype(AMREX_D_TERM(i, | (j<<1), | (k<<2)))
        {}
    // dtor, copy-ctor, copy-op=, move-ctor, and move-op= are compiler generated.

    //! Set IndexType to be NODE based in direction dir.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void set (int dir) noexcept { itype |= mask(dir); }
    //! Set IndexType to be CELL based in direction dir.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void unset (int dir) noexcept { itype &= ~mask(dir); }
    //! True if IndexType is NODE based in direction dir.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool test (int dir) const noexcept { return (itype & mask(dir)) != 0; }
    //! Set NODE based in all directions.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void setall () noexcept { itype = (1 << AMREX_SPACEDIM) - 1; }
    //! Set CELL based in all directions.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void clear () noexcept { itype = 0; }
    //! True if this IndexType is NODE based in any direction.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool any () const noexcept { return itype != 0; }
    //! True if IndexType is valid.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool ok () const noexcept { return itype < (1 << AMREX_SPACEDIM); }
    //! Change from CELL to NODE or NODE to CELL in direction dir.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void flip (int i) noexcept { itype ^= mask(i); }
    //! True if IndexTypes are identical.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool operator== (const IndexType& t) const noexcept { return t.itype == itype; }
    //! True if IndexTypes are not identical.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool operator!= (const IndexType& t) const noexcept { return t.itype != itype; }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool operator< (const IndexType& t) const noexcept { return itype < t.itype; }
    //! True if the IndexType is CELL based in all directions.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool cellCentered () const noexcept { return itype == 0; }
    //! True if the IndexType is CELL based in dir-direction.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool cellCentered (int dir) const noexcept { return (itype & mask(dir)) == 0; }
    //! True if the IndexType is NODE based in all directions.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool nodeCentered () const noexcept { return itype == (1<<AMREX_SPACEDIM)-1; }
    //! True if the IndexType is NODE based in dir-direction.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool nodeCentered (int dir) const noexcept { return (itype & mask(dir)) != 0; }
    //! Set IndexType to CellIndex type t in direction dir.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void setType (int dir, CellIndex t) noexcept { t == CELL ? unset(dir) : set(dir); }
    //! Returns  the CellIndex in direction dir.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    CellIndex ixType (int dir) const noexcept { return (CellIndex) ((itype & (1<<dir)) >> dir); }
    //! Return an integer representing the IndexType in direction dir.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    int operator[] (int dir) const noexcept { return test(dir); }
    //! Fill an IntVect of size AMREX_SPACEDIM with IndexTypes.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect ixType () const noexcept { return IntVect(AMREX_D_DECL(itype&1, (itype>>1)&1, (itype>>2)&1)); }
    /**
    * \brief This static member function returns an IndexType object of value
    * IndexType::CELL.  It is provided as a convenience to our users
    * when defining a Box all of whose faces should be of type
    * IndexType::CELL.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static constexpr IndexType TheCellType () noexcept {
        return IndexType(AMREX_D_DECL(IndexType::CELL,
                                      IndexType::CELL,
                                      IndexType::CELL));
    }
    /**
    * \brief This static member function returns an IndexType object of value
    * IndexType::NODE.  It is provided as a convenience to our users
    * when defining a Box all of whose faces should be of type
    * IndexType::NODE.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static constexpr IndexType TheNodeType () noexcept {
        return IndexType(AMREX_D_DECL(IndexType::NODE,
                                      IndexType::NODE,
                                      IndexType::NODE));
    }

private:
    //! Returns 1<<k.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static int mask (int k) noexcept { return 1<<k; }
    //! An integer holding the CellIndex in bits 0 - AMREX_SPACEDIM-1.
    unsigned int itype;
};

//! Write an IndexType to an ostream in ASCII.
std::ostream& operator<< (std::ostream& os, const IndexType& itype);
//! Read an IndexType from an istream.
std::istream& operator>> (std::istream& is, IndexType& itype);

}

#endif /*BL_INDEXTYPE_H*/
