
#ifndef AMREX_BNDRYDATA_H_
#define AMREX_BNDRYDATA_H_
#include <AMReX_Config.H>

#include <AMReX_BndryRegister.H>
#include <AMReX_BoundCond.H>
#include <AMReX_MultiMask.H>

#include <memory>
#include <map>

namespace amrex {

/**
* \brief A BndryData stores and manipulates boundary data information on
* each side of each box in a BoxArray.
*
*    A BndryData contains a BndryRegister about each side of each grid in
*    a Boxarray.  These data are used to store information along the
*    outer layer of each grid (at the same level of coarseness), such
*    as the value of boundary conditions there.  Also, for each
*    boundary, this class maintains a BoundCond identifier and a
*    location.  This "location" often specifies where, in physical
*    space, a value, which is stored in the boundary registers, will
*    be applied by a solver (although interpretation of this information
*    is application-specific).
*
*    In addition to boundary condition values, types and locations, and
*    BndryDate object maintains a mask for each boundary value.  Given
*    a Geometry describing the physical "domain" of application, a BndryData
*    object fills the mask arrays with a mask value to indicate if the
*    node is outside the domain (outside_domain), or, if not, whether it
*    lays within the valid region of an adjacent grid (covered), or not
*    (not_covered).  This mask data is created upon non-default class
*    instantiation.
*/
template <typename MF>
class BndryDataT
    : public BndryRegisterT<MF>
{
public:
    //! mask values enumeration
    enum MaskVal { covered = 0, not_covered = 1, outside_domain = 2, NumMaskVals = 3 };

    using value_type = typename MF::value_type;

    //! Default constructor
    BndryDataT() noexcept = default;

    /**
    * \brief constructor specifying number of components and box of physical
    * domain (cell-centered)
    */
    BndryDataT (const BoxArray& grids,
                const DistributionMapping& dmap,
                int             ncomp,
                const Geometry& geom);

    //! destructor
    ~BndryDataT () = default;

    BndryDataT (const BndryDataT<MF>& src) = delete;
    BndryDataT (BndryDataT<MF>&& rhs) = delete;
    BndryDataT<MF>& operator= (const BndryDataT<MF>& src) = delete;
    BndryDataT<MF>& operator= (BndryDataT<MF>&& rhs) = delete;

    //! allocate bndry fabs along given face
    void define (const BoxArray& grids,
                 const DistributionMapping& dmap,
                 int             ncomp,
                 const Geometry& geom);
    //
    const MultiMask& bndryMasks (Orientation face) const noexcept { return masks[face]; }

    //! Return FabSet on given face.
    const FabSetT<MF>& bndryValues (Orientation face) const noexcept {
        return this->bndry[face];
    }

    //! Some useful typedefs
    using RealTuple = Array<value_type,2*AMREX_SPACEDIM>;

    /**
    * \brief Return boundary location on given face for grids we own.
    * It's an error if we don't own that grid.
    * RealTuple is indexed with Orientation.
    */
    const RealTuple& bndryLocs (int igrid) const noexcept;
    const RealTuple& bndryLocs (const MFIter& mfi) const noexcept;

    /**
    * \brief Return boundary type specifier on given face for grids we own.
    * It's an error if we don't own that grid.
    */
    const Vector< Vector<BoundCond> >& bndryConds (int igrid) const noexcept;
    const Vector< Vector<BoundCond> >& bndryConds (const MFIter& mfi) const noexcept;

    //! return number of components for which this object is intended
    int nComp () const noexcept { return m_ncomp; }

    //! return domain used to define masks
    const Box& getDomain () const noexcept { return geom.Domain(); }

    //! return geometry used to define masks
    const Geometry& getGeom () const noexcept { return geom; }

    //! set values of boundary Fab for given orientation on nth grid
    void setValue (Orientation face, int n, Real val) noexcept;

    //! set boundary type specifier for given orientation on nth grid
    void setBoundCond (Orientation     face,
                       int              n,
                       int              comp,
                       const BoundCond& bcn) noexcept;

    void setBoundCond (Orientation     face,
                       const MFIter&    mfi,
                       int              comp,
                       const BoundCond& bcn) noexcept;

    //! set boundary location for given orientation on nth grid
    void setBoundLoc (Orientation face,
                      int         n,
                      value_type val) noexcept;

    void setBoundLoc (Orientation face,
                      const MFIter& mfi,
                      value_type val) noexcept;

    //! implement public access to const BndryRegister::operator[]
    const FabSetT<MF>& operator[] (Orientation face) const noexcept { return BndryRegisterT<MF>::bndry[face]; }

    //! implement public access to BndryRegister::operator[]
    FabSetT<MF>& operator[] (Orientation face) noexcept { return BndryRegisterT<MF>::bndry[face]; }

protected:
    /**
    * \brief Map of boundary condition type specifiers.
    * The outer Array dimension is over Orientation.
    */
    LayoutData<Vector< Vector<BoundCond> > > bcond;

    LayoutData<RealTuple> bcloc;
    //! Boundary condition mask
    Vector<MultiMask> masks;
    //! Domain used for mask definitions.
    Geometry geom;
    int      m_ncomp{-1};
    bool     m_defined{false};

private:
    // Mask info required for this many cells past grid edge (here,
    //  e.g. ref=4, crse stencil width=3, and let stencil slide 2 past grid
    //  edge)
    static constexpr int NTangHalfWidth = 5; // ref_ratio + 1, so won't work if ref_ratio > 4
};

template <typename MF>
BndryDataT<MF>::BndryDataT (const BoxArray& _grids,
                            const DistributionMapping& _dmap,
                            int             _ncomp,
                            const Geometry& _geom)
    :
    geom(_geom),
    m_ncomp(_ncomp)
{
    define(_grids,_dmap,_ncomp,_geom);
}

template <typename MF>
void
BndryDataT<MF>::setBoundCond (Orientation     _face,
                              int              _n,
                              int              _comp,
                              const BoundCond& _bcn) noexcept
{
    bcond[_n][_face][_comp] = _bcn;
}

template <typename MF>
void
BndryDataT<MF>::setBoundCond (Orientation     _face,
                              const MFIter&   mfi,
                              int              _comp,
                              const BoundCond& _bcn) noexcept
{
    bcond[mfi][_face][_comp] = _bcn;
}

template <typename MF>
void
BndryDataT<MF>::setBoundLoc (Orientation _face,
                             int         _n,
                             value_type _val) noexcept
{
    bcloc[_n][_face] = _val;
}

template <typename MF>
void
BndryDataT<MF>::setBoundLoc (Orientation _face,
                             const MFIter& mfi,
                             value_type _val) noexcept
{
    bcloc[mfi][_face] = _val;
}

template <typename MF>
const Vector< Vector<BoundCond> >&
BndryDataT<MF>::bndryConds (int igrid) const noexcept
{
    return bcond[igrid];
}

template <typename MF>
const Vector< Vector<BoundCond> >&
BndryDataT<MF>::bndryConds (const MFIter& mfi) const noexcept
{
    return bcond[mfi];
}

template <typename MF>
const typename BndryDataT<MF>::RealTuple&
BndryDataT<MF>::bndryLocs (int igrid) const noexcept
{
    return bcloc[igrid];
}

template <typename MF>
const typename BndryDataT<MF>::RealTuple&
BndryDataT<MF>::bndryLocs (const MFIter& mfi) const noexcept
{
    return bcloc[mfi];
}

template <typename MF>
void
BndryDataT<MF>::define (const BoxArray& _grids,
                        const DistributionMapping& _dmap,
                        int             _ncomp,
                        const Geometry& _geom)
{
    BL_PROFILE("BndryData::define()");

    if (m_defined) {
        if (_grids == this->boxes() && m_ncomp == _ncomp && _geom.Domain() == geom.Domain()) {
            // We want to allow reuse of BndryDataT objects that were define()d exactly as a previous call.
            return;
        }
        // Otherwise we'll just abort.  We could make this work but it's just as easy to start with a fresh Bndrydata object.
        amrex::Abort("BndryDataT<MF>::define(): object already built");
    }
    geom    = _geom;
    m_ncomp = _ncomp;

    BndryRegisterT<MF>::setBoxes(_grids);

    masks.clear();
    masks.resize(2*AMREX_SPACEDIM);

    for (OrientationIter fi; fi; ++fi) {
        Orientation face = fi();
        BndryRegisterT<MF>::define(face,IndexType::TheCellType(),0,1,1,_ncomp,_dmap);
        masks[face].define(_grids, _dmap, geom, face, 0, 2, NTangHalfWidth, 1, true);
    }
    //
    // Define "bcond" and "bcloc".
    //
    // We note that all orientations of the FabSets have the same distribution.
    // We'll use the low 0 side as the model.
    //
    //

    bcloc.define(_grids, _dmap);
    bcond.define(_grids, _dmap);

    for (FabSetIter bfsi(this->bndry[Orientation(0,Orientation::low)]);
         bfsi.isValid(); ++bfsi) {
        Vector< Vector<BoundCond> >& abc = bcond[bfsi];
        abc.resize(2*AMREX_SPACEDIM);
        for (OrientationIter fi; fi; ++fi) {
            abc[fi()].resize(_ncomp);
        }
    }

    m_defined = true;
}

template <typename MF>
void
BndryDataT<MF>::setValue (Orientation face, int n, Real val) noexcept
{
    auto& fab = this->bndry[face][n];
    auto arr = fab.array();
    const Box& bx = fab.box();
    const int ncomp = m_ncomp;
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, m,
    {
        arr(i,j,k,m) = val;
    });
}

using BndryData = BndryDataT<MultiFab>;
using fBndryData = BndryDataT<fMultiFab>;

}

#endif /*_BNDRYDATA_H_*/
