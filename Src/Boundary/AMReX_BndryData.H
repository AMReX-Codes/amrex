
#ifndef AMREX_BNDRYDATA_H_
#define AMREX_BNDRYDATA_H_

#include <AMReX_BLassert.H>
#include <AMReX_BoxArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Orientation.H>
#include <AMReX_BndryRegister.H>
#include <AMReX_MultiMask.H>
#include <AMReX_BoundCond.H>
#include <AMReX_Geometry.H>
#include <AMReX_Array.H>
#include <AMReX_LayoutData.H>

#include <memory>
#include <map>

namespace amrex {

//
//  A BndryData stores and manipulates boundary data
//  information on each side of each box in a BoxArray.
//    

/**
* \brief A BndryData contains a BndryRegister about each side of each grid in
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

class BndryData
    :
    public BndryRegister
{
public:
    //! mask values enumeration
    enum MaskVal { covered = 0, not_covered = 1, outside_domain = 2, NumMaskVals = 3 };

    //! Default constructor
    BndryData() noexcept;

    /**
    * \brief constructor specifying number of components and box of physical
    * domain (cell-centered)
    */
    BndryData (const BoxArray& grids,
	       const DistributionMapping& dmap,
               int             ncomp,
               const Geometry& geom);

    //! destructor
    virtual ~BndryData ();

    //! copy constructor
    BndryData (const BndryData& src);

    //! copy assignment operator
    BndryData& operator= (const BndryData& src);

    BndryData (BndryData&& rhs) = delete;
    BndryData& operator= (BndryData&& rhs) = delete;

    //! alocate bndry fabs along given face
    void define (const BoxArray& grids,
		 const DistributionMapping& dmap,
                 int             ncomp,
                 const Geometry& geom);
    //
    const MultiMask& bndryMasks (Orientation face) const noexcept { return masks[face]; }

    //! Return FabSet on given face.
    const FabSet& bndryValues (Orientation face) const noexcept { return bndry[face]; }

    //! Some useful typedefs
    using RealTuple = Array<Real,2*AMREX_SPACEDIM>;

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
                      Real        val) noexcept;

    void setBoundLoc (Orientation face,
                      const MFIter& mfi,
                      Real        val) noexcept;

    //! implement public access to const BndryRegister::operator[]
    const FabSet& operator[] (Orientation face) const noexcept { return BndryRegister::bndry[face]; }

    //! implement public access to BndryRegister::operator[]
    FabSet& operator[] (Orientation face) noexcept { return BndryRegister::bndry[face]; }

protected:
    //! Helper function for copy constructor and assigment operator.
    void init (const BndryData& src);

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
    int      m_ncomp;
    bool     m_defined;

private:
    static int NTangHalfWidth;
};

}

#endif /*_BNDRYDATA_H_*/

