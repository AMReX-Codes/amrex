
#ifndef _BNDRYDATA_H_
#define _BNDRYDATA_H_

#include <BLassert.H>
#include <BoxArray.H>
#include <FArrayBox.H>
#include <Orientation.H>
#include <BndryRegister.H>
#include <Mask.H>
#include <BoundCond.H>
#include <Geometry.H>
#include <Tuple.H>

#include <map>

//
//  A BndryData stores and manipulates boundary data
//  information on each side of each box in a BoxArray.
//    

/*
    A BndryData contains a BndryRegister about each side of each grid in
    a Boxarray.  These data are used to store information along the
    outer layer of each grid (at the same level of coarseness), such
    as the value of boundary conditions there.  Also, for each
    boundary, this class maintains a BoundCond identifier and a
    location.  This "location" often specifies where, in physical
    space, a value, which is stored in the boundary registers, will
    be applied by a solver (although interpretation of this information
    is application-specific).

    In addition to boundary condition values, types and locations, and
    BndryDate object maintains a mask for each boundary value.  Given
    a Geometry describing the physical "domain" of application, a BndryData
    object fills the mask arrays with a mask value to indicate if the
    node is outside the domain (outside_domain), or, if not, whether it
    lays within the valid region of an adjacent grid (covered), or not
    (not_covered).  This mask data is created upon non-default class
    instantiation.
*/

class BndryData
    :
    private BndryRegister
{
public:
    //
    // return the array of Boxes
    //
    using BndryRegister::boxes;
    //
    // return the number of Boxes
    //
    using BndryRegister::size;
    //
    // mask values enumeration
    //
    enum MaskVal { covered = 0, not_covered = 1, outside_domain = 2, NumMaskVals = 3 };
    //
    // Default constructor
    //
    BndryData();
    //
    // constructor specifying number of components and box of physical
    // domain (cell-centered)
    //
    BndryData (const BoxArray& grids,
               int             ncomp,
               const Geometry& geom,
	       ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());
    //
    // destructor
    //
    virtual ~BndryData ();
    //
    // copy constructor
    //
    BndryData (const BndryData& src);
    //
    // copy assignment operator
    //
    BndryData& operator= (const BndryData& src);
    //
    // alocate bndry fabs along given face
    //
    void define (const BoxArray& grids,
                 int             ncomp,
                 const Geometry& geom,
		 ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());
    //
    // write to output stream
    //
    friend std::ostream& operator<< (std::ostream& os, const BndryData& bd);
    //
    // Return FabSet on given face.
    //
    const FabSet& bndryValues (Orientation face) const { return bndry[face]; }
    //
    // Some useful typedefs
    //
    typedef Tuple<Real,2*BL_SPACEDIM> RealTuple;

    typedef Tuple<Mask*,2*BL_SPACEDIM> MaskTuple;
    //
    // Return boundary location on given face for grids we own.
    // It's an error if we don't own that grid.
    // RealTuple is indexed with Orientation.
    //
    const RealTuple& bndryLocs (int igrid) const;
    //
    // Return boundary type specifier on given face for grids we own.
    // It's an error if we don't own that grid.
    //
    const Array< Array<BoundCond> >& bndryConds (int igrid) const;
    //
    // Return boundary masks for grids we own.
    // It's an error if we don't own that grid.
    //
    const MaskTuple& bndryMasks (int igrid) const;
    //
    // return number of components for which this object is intended
    //
    int nComp () const { return m_ncomp; }
    //
    using BndryRegister::color;
    //
    // return domain used to define masks
    //
    const Box& getDomain () const { return geom.Domain(); }
    //
    // return geometry used to define masks
    //
    const Geometry& getGeom () const { return geom; }
    //
    // set values of boundary Fab for given orientation on nth grid
    //
    void setValue (Orientation face, int n, Real val) { bndry[face][n].setVal(val); }
    //
    // set mask values for given orientation on nth grid
    //
    void setMaskValue (Orientation face,
                       int         n,
                       int         val);
    //
    // set boundary type specifier for given orientation on nth grid
    //
    void setBoundCond (Orientation     face,
                       int              n,
                       int              comp,
                       const BoundCond& bcn);
    //
    // set boundary location for given orientation on nth grid
    //
    void setBoundLoc (Orientation face,
                      int         n,
                      Real        val);
    //
    // implement public access to const BndryRegister::operator[]
    //
    const FabSet& operator[] (Orientation face) const { return BndryRegister::bndry[face]; }
    //
    // implement public access to BndryRegister::operator[]
    //
    FabSet& operator[] (Orientation face) { return BndryRegister::bndry[face]; }
    //
    // write the boundary data object to an ostream
    //
    void writeOn (std::ostream& os) const;
    //
    // read the boundary data object from an istream
    //
    void readFrom (std::istream& is);

protected:
    //
    // Helper function for copy constructor and assigment operator.
    //
    void init (const BndryData& src);
    //
    // protect BndryRegister grids.
    //
    using BndryRegister::grids;
    //
    // protect BndryRegister FabSets.
    //
    using BndryRegister::bndry;
    //
    // Map of boundary condition type specifiers.
    // The outer Array dimension is over Orientation.
    //
    std::map< int, Array< Array<BoundCond> > > bcond;

    std::map<int,RealTuple> bcloc;
    //
    // Map boundary condition mask arrays (FAB-like).
    //
    std::map<int,MaskTuple> masks;
    //
    // Domain used for mask definitions.
    //
    Geometry geom;
    int      m_ncomp;
    bool     m_defined;

private:
    //
    // Free memory taken by masks for this BndryData object.
    //
    void clear_masks ();

    static int NTangHalfWidth;
};

#endif /*_BNDRYDATA_H_*/

