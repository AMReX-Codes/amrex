//BL_COPYRIGHT_NOTICE

#ifndef _BNDRYDATA_H_
#define _BNDRYDATA_H_

//
// $Id: BndryData.H,v 1.4 1998-07-07 17:25:12 lijewski Exp $
//

#include <Assert.H>
#include <BoxArray.H>
#include <FArrayBox.H>
#include <Orientation.H>
#include <BndryRegister.H>
#include <Mask.H>
#include <BoundCond.H>
#include <Geometry.H>

//@Man:
/*@Memo:
  A ProxyGeometry is a class constructed solely to prevent automatic
  type conversion with the Geometry::Geometry(const Box&) constructor.
  This class can easily make a Geometry object, but does not have a
  ProxyGeometry::ProxyGeometry(const Box&) constructor.
*/    
/*@Doc:
        A ProxyGeometry object is, for most purposes, merely a wrapper to the
    Geometry class.  The Geometry class contains a single-argument
    constructor taking a const Box &, and therefore allows implicit
    type conversion.  However, that constructor requires parm-parsed
    input for the coordinate system and prob lo and hi.  This has led to
    some some mysterious error messages, by correctly compiled code
    passing in a Box where a Geometry was called for.  This little class
    is an attempt to circumvent this little funny business.
*/

class ProxyGeometry
{

public:
    //
    //@ManDoc: The default constructor.
    //
    ProxyGeometry (const Geometry& g)
        :
        m_g(g) {}
    //
    //@ManDoc: Type conversion operator.
    //
    operator Geometry () const { return m_g; }

protected:
    //
    //@ManDoc: The underlying Geometry
    //
    const Geometry& m_g;
};

//@Man:
/*@Memo:
  A BndryData stores and manipulates boundary data
  information on each side of each box in a BoxArray.
*/    
/*@Doc:
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
    //@ManDoc: return the array of Boxes
    //
    BndryRegister::boxes;
    //
    //@ManDoc: return the number of Boxes
    //
    BndryRegister::length;
    //
    //@ManDoc: mask values enumeration
    //
    enum MaskVal { covered = 0, not_covered = 1, outside_domain = 2, NumMaskVals = 3 };
    //
    //@ManDoc: default constructor
    //
    BndryData() : BndryRegister () {}

    /*@ManDoc: constructor specifying number of components and box of physical
        domain (cell-centered)
    */
    BndryData (const BoxArray&     _grids,
               int                  _ncomp,
               const ProxyGeometry& geom);
    //
    //@ManDoc: destructor
    //
    virtual ~BndryData ();
    //
    //@ManDoc: copy constructor
    //
    BndryData (const BndryData& src);
    //
    //@ManDoc: copy operator
    //
    BndryData& operator= (const BndryData& src);
    //
    //@ManDoc: alocate bndry fabs along given face
    //
    void define (const BoxArray&      _grids,
                 int                  _ncomp,
                 const ProxyGeometry& geom);
    //
    //@ManDoc: write to output stream
    //
    friend ostream& operator<< (ostream&         os,
                                const BndryData& bd);
    //
    //@ManDoc: return FabSet on given face
    //
    const FabSet& bndryValues (const Orientation& _face) const;
    //
    //@ManDoc: return boundary location on given face
    //
    const Array<Real>& bndryLocs (const Orientation& _face) const;
    //
    //@ManDoc: return boundary type specifyer on given face
    //
    const Array<Array<BoundCond> >& bndryConds (const Orientation& _face) const;
    //
    //@ManDoc: return boundary mask on given face
    //
    const PArray<Mask>& bndryMasks (const Orientation& _face) const;
    //
    //@ManDoc: return domain used to define masks
    //
    const Box& getDomain () const;
    //
    //@ManDoc: return geometry used to define masks
    //
    const Geometry& getGeom () const;
    //
    //@ManDoc: set values of boundary Fab for given orientation on nth grid
    //
    void setValue (const Orientation& _face,
                   int                _n,
                   Real               _val);
    //
    //@ManDoc: set mask values for given orientation on nth grid
    //
    void setMaskValue (const Orientation& _face,
                       int                _n,
                       int                _val);
    //
    //@ManDoc: set boundary type specifier for given orientation on nth grid
    //
    void setBoundCond (const Orientation& _face,
                       int                _n,
                       int                _comp,
                       const BoundCond&   _bcn);
    //
    //@ManDoc: set boundary location for given orientation on nth grid
    //
    void setBoundLoc (const Orientation& _face,
                      int                _n,
                      Real               _val);
    //
    //@ManDoc: implement public access to const BndryRegister::operator[]
    //
    const FabSet& operator[] (const Orientation& _face) const;
    //
    //@ManDoc: implement public access to BndryRegister::operator[]
    //
    FabSet& operator[] (const Orientation& _face);

protected:
    
    /*@ManDoc: Array (on orientation)(on grid)(on component) of boundary
         condition type specifiers.
    */
    Array< Array< BoundCond > > bcond[2*BL_SPACEDIM];
    //
    //@ManDoc: Array (on orientation) of boundary condition locations.
    //
    Array< Real > bcloc[2*BL_SPACEDIM];
    //
    //@ManDoc: Array (on orientation) of boundary condition mask arrays (FAB-like)
    //
    PArray< Mask > masks[2*BL_SPACEDIM];
    //
    //@ManDoc: Domain used for mask defintions.
    //
    Geometry geom;
    //
    //@ManDoc: protect BndryRegister grids.
    //
    BndryRegister::grids;
    //
    //@ManDoc: protect BndryRegister FabSets.
    //
    BndryRegister::bndry;

private:
    //
    // Free memory taken by masks for this BndryData object.
    //
    void clear_masks ();
};

inline
const FabSet&
BndryData::bndryValues (const Orientation& _face) const
{
    return bndry[_face];
}

inline
const Array<Real>&
BndryData::bndryLocs (const Orientation& _face) const
{
    return bcloc[_face];
}

inline
const Array<Array<BoundCond> >&
BndryData::bndryConds (const Orientation& _face) const
{
    return bcond[_face];
}

inline
const PArray<Mask>&
BndryData::bndryMasks (const Orientation& _face) const
{
    return masks[_face];
}

inline
const Box&
BndryData::getDomain () const
{
    return geom.Domain();
}

inline
const Geometry&
BndryData::getGeom () const
{
    return geom;
}

inline
void
BndryData::setValue (const Orientation& _face,
                     int                _n,
                     Real               _val)
{
    bndry[_face][_n].setVal(_val);
}

inline
void
BndryData::setMaskValue (const Orientation& _face,
                         int                _n,
                         int                _val)
{
    masks[_face][_n].setVal(_val);
}

inline
void
BndryData::setBoundCond (const Orientation& _face,
                         int                _n,
                         int                _comp,
                         const BoundCond&   _bcn)
{
    bcond[_face][_n][_comp] = _bcn;
}

inline
void
BndryData::setBoundLoc (const Orientation& _face,
                        int                _n,
                        Real               _val)
{
    bcloc[_face][_n] = _val;
}

inline
const FabSet&
BndryData::operator[] (const Orientation& _face) const
{
    return BndryRegister::bndry[_face];
}

inline
FabSet&
BndryData::operator[] (const Orientation& _face)
{
    return BndryRegister::bndry[_face];
}

#endif /*_BNDRYDATA_H_*/

