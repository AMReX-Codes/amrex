//BL_COPYRIGHT_NOTICE

#ifndef _COORDSYS_H_
#define _COORDSYS_H_ 

//
// $Id: CoordSys.H,v 1.4 1998-05-06 22:36:15 lijewski Exp $
//

#include <REAL.H>
#include <Array.H>
#include <Box.H>

//
// Forward declaration.
//
class FArrayBox;

//
//@Man:
//@Memo: Coordinate System
/*@Doc:

  Routines for mapping between physical coordinate system and index space.
*/

class CoordSys
{
public:
#if (BL_SPACEDIM == 2)    
    enum CoordType { undef = -1, cartesian = 0, RZ = 1 };
#endif
#if (BL_SPACEDIM == 3)
    enum CoordType { undef = -1, cartesian = 0 };
#endif    
    //
    //@ManDoc: Nice ASCII output.
    //
    friend ostream& operator<< (ostream&, const CoordSys& );
    //
    //@ManDoc: Nice ASCII input.
    //
    friend istream& operator>> (istream&, CoordSys& );
    //
    //@ManDoc: Default constructor to undefined state.
    //
    CoordSys ();
    //
    //@ManDoc: Construct, specify a cell size, and low end of index range.
    //
    CoordSys (const Real* cell_dx);
    //
    //@ManDoc: The destructor.
    //
    ~CoordSys ();
    //
    //@ManDoc: Initialize after construction.
    //
    void define (const Real* cell_dx);
    //
    //@ManDoc: Set the CoordType.
    //
    static void SetCoord (CoordType coord);
    //
    //@ManDoc: Returns the CoordType.
    //
    static CoordType Coord ();
    //
    //@ManDoc: Is CoordType == RZ?
    //
    static bool IsRZ ();
    //
    //@ManDoc: Is CoordType == cartesion?
    //
    static bool IsCartesian ();
    //
    //@ManDoc: Sets the offset for each coordinate direction.
    //
    static void SetOffset (const Real* x_lo);
    //
    //@ManDoc: Returns the offset.
    //
    static const Real* Offset ();
    //
    //@ManDoc: Returns the offset for the specified coordinate direction.
    //
    static Real Offset (int dir);
    //
    //@ManDoc: Returns the cellsize for each coordinate direction.
    //
    const Real* CellSize () const;
    //
    //@ManDoc: Returns the cellsize for the specified coordinate direction.
    //
    Real CellSize (int dir) const;
    //
    //@ManDoc: Returns location of cell center in specified direction.
    //
    Real CellCenter (int point, int dir) const;
    //
    //@ManDoc: Return location of cell center.
    //
    void CellCenter (const IntVect& point, Array<Real>& loc) const;
    //
    //@ManDoc: Return location of cell center.
    //
    void CellCenter (const IntVect& point, Real* loc) const;
    //
    //@ManDoc: Returns location of lo edge in specified direction.
    //
    Real LoEdge (int point, int dir) const;
    //
    //@ManDoc: Equivalent to LoEdge(point[dir], dir).
    //
    Real LoEdge (const IntVect& point, int dir) const;
    //
    //@ManDoc: Returns location of hi edge in specified direction.
    //
    Real HiEdge (int point, int dir) const;
    //
    //@ManDoc: Equivalent to HiEdge(point[dir], dir).
    //
    Real HiEdge (const IntVect& point, int dir) const;
    //
    //@ManDoc: Sets location of lo face into `loc'.
    //
    void LoFace (const IntVect& point, int dir, Array<Real>& loc) const;
    //
    //@ManDoc: Sets location of lo face into `loc'.
    //
    void LoFace (const IntVect& point, int dir, Real* loc) const;
    //
    //@ManDoc: Sets location of hi face into `loc'.
    //
    void HiFace (const IntVect& point, int dir, Array<Real>& loc) const;
    //
    //@ManDoc: Sets location of hi face into `loc'.
    //
    void HiFace (const IntVect& point, int dir, Real* loc) const;
    //
    //@ManDoc: Return location of lower left hand corner.
    //
    void LoNode (const IntVect& point, Array<Real>& loc) const;
    //
    //@ManDoc: Return location of lower left hand corner.
    //
    void LoNode (const IntVect& point, Real* loc) const;
    //
    //@ManDoc: Return location of upper right hand corner.
    //
    void HiNode (const IntVect& point, Array<Real>& loc) const;
    //
    //@ManDoc: Return location of upper right hand corner.
    //
    void HiNode (const IntVect& point, Real* loc) const;

    /*@ManDoc: Returns cell centered index of cell containing point.
               This may return undesired results if point
               is on a cell boundary.
    */ 
    IntVect CellIndex (const Real* point) const;

    /*@ManDoc: Returns node centered index of lower left hand corner of
               cell containing this point.
    */
    IntVect LowerIndex (const Real* point) const;

    /*@ManDoc: Returns node centered index of upper right hand corner of
               cell containing this point.
    */
    IntVect UpperIndex (const Real* point) const;

    /*@ManDoc: Compute cell volumes in given region and place them into
               resize()d input FAB.
    */
    void GetVolume (FArrayBox& vol, const Box& region) const;

    /*@ManDoc: Compute cell volumes in given region and place them into
               new()d FAB.  It is the user's respoinsibility to delete
               the returned FAB.
    */
    FArrayBox* GetVolume (const Box& region) const;
    
#if (BL_SPACEDIM == 2)
    /*@ManDoc: Compute d(log(A))/dr at cell centers in given region and
               return the results in the resize()d input FAB.
    */
    void GetDLogA (FArrayBox& dloga, const Box& region, int dir) const;

    /*@ManDoc: Compute d(log(A))/dr at cell centers in given region and
               return the results in a new()d FAB.  It is the user's
               responsibility to delete the returned FAB.
    */
    FArrayBox* GetDLogA (const Box& region, int dir) const;
#endif
    //
    //@ManDoc: Return the volume of the specified cell.
    //
    Real Volume (const IntVect& point) const;
    //
    //@ManDoc: Return the volume of the specified cell.
    //
    Real Volume (const Real xlo[BL_SPACEDIM],
                 const Real xhi[BL_SPACEDIM]) const;

    /*@ManDoc: Compute area of cell faces in given region and given
               index direction and return the result in resize()d input FAB.
    */
    void GetFaceArea (FArrayBox& area, const Box& region, int dir) const;

    /*@ManDoc: Compute area of cell faces in given region and given
               index direction and return the result in new()d FAB.
               It is the user's responsibility to delete the FAB.
    */
    FArrayBox* GetFaceArea (const Box& region, int dir) const;
    //
    //@ManDoc: Returns lo face area of given cell in direction dir.
    //
    Real AreaLo (const IntVect& point, int dir) const;
    //
    //@ManDoc: Returns hi face area of given cell in direction dir.
    //
    Real AreaHi (const IntVect& point, int dir) const;

    /*@ManDoc: Return array of physical locations of cell edges
               in the resize()d input array.
    */
    void GetEdgeLoc (Array<Real>& loc, const Box& region, int dir) const;

    /*@ManDoc: Return array of physical locations of cell centers
               in the resize()d input array.
    */
    void GetCellLoc (Array<Real>& loc, const Box& region, int dir) const;

    /*@ManDoc: Return array of volume coordinates at cell edges
               for region in given direction.
    */
    void GetEdgeVolCoord (Array<Real>& loc, const Box& region, int dir) const;

    /*@ManDoc: Return array of volume coordinates at cell centers
               for region in given direction.
    */
    void GetCellVolCoord (Array<Real>& loc, const Box& region, int dir) const;

protected:
    //
    // Static data.
    //
    static CoordType c_sys;
    static Real      offset[BL_SPACEDIM];
    //
    // Non-static data.
    //  
    Real dx[BL_SPACEDIM];
    bool ok;
};

//
// Inlines.
//

inline CoordSys::~CoordSys () {}

inline
CoordSys::CoordType
CoordSys::Coord ()
{
    return c_sys;
}

inline
const Real*
CoordSys::Offset ()
{
    return offset;
}

inline
Real
CoordSys::Offset (int dir)
{
    return offset[dir];
}

inline
const Real*
CoordSys::CellSize () const
{
    assert(ok);
    return dx;
}

inline
Real
CoordSys::CellSize (int dir) const
{
    assert(ok);
    return dx[dir];
}

inline
void
CoordSys::SetCoord (CoordType coord)
{
    c_sys = coord;
}

inline
void
CoordSys::SetOffset (const Real* x_lo)
{
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        offset[k] = x_lo[k];
    }
}

inline
bool
CoordSys::IsRZ ()
{
    assert(c_sys != undef);
#if (BL_SPACEDIM == 2)
    return (c_sys == RZ);
#endif    
#if (BL_SPACEDIM == 3)
    return false;
#endif    
}

inline
bool
CoordSys::IsCartesian ()
{
    assert(c_sys != undef);
    return (c_sys == cartesian);
}

inline
CoordSys::CoordSys ()
{
    ok = false;
}

inline
void
CoordSys::define (const Real* cell_dx)
{
    assert(c_sys != undef);
    ok = true;
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        dx[k] = cell_dx[k];
    }
}

inline
CoordSys::CoordSys (const Real* cell_dx)
{
    define(cell_dx);
}

inline
void
CoordSys::CellCenter (const IntVect& point,
                      Real*          loc) const
{
    assert(ok);
    assert(loc != 0);
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        loc[k] = offset[k] + dx[k]*(0.5+ (Real)point[k]);
    }
}

inline
void
CoordSys::CellCenter (const IntVect& point,
                      Array<Real>&   loc) const
{
    assert(ok);
    loc.resize(BL_SPACEDIM);
    CellCenter(point, loc.dataPtr());
}

inline
Real
CoordSys::CellCenter (int point,
                      int dir) const
{
    assert(ok);
    return offset[dir] + dx[dir]*(0.5+ (Real)point);
}

inline
Real
CoordSys::LoEdge (int point,
                  int dir) const
{
    assert(ok);
    return offset[dir] + dx[dir]*point;
}

inline
Real
CoordSys::LoEdge (const IntVect& point,
                  int            dir) const
{
    assert(ok);
    return offset[dir] + dx[dir]*point[dir];
}

inline
Real
CoordSys::HiEdge (int point,
                  int dir) const
{
    assert(ok);
    return offset[dir] + dx[dir]*(point + 1);
}

inline
Real
CoordSys::HiEdge (const IntVect& point,
                  int            dir) const
{
    assert(ok);
    return offset[dir] + dx[dir]*(point[dir] + 1);
}

inline
void
CoordSys::LoFace (const IntVect& point,
                  int            dir,
                  Real*          loc) const
{
    assert(ok);
    assert(loc != 0);
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        Real off = (k == dir) ? 0.0 : 0.5;
        loc[k] = offset[k] + dx[k]*(off + (Real)point[k]);
    }
}

inline
void
CoordSys::LoFace (const IntVect& point,
                  int            dir,
                  Array<Real>&   loc) const
{
    loc.resize(BL_SPACEDIM);
    LoFace(point,dir, loc.dataPtr());
}

inline
void
CoordSys::HiFace (const IntVect& point,
                  int            dir,
                  Real*          loc) const
{
    assert(ok);
    assert(loc != 0);
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        Real off = (k == dir) ? 1.0 : 0.5;
        loc[k] = offset[k] + dx[k]*(off + (Real)point[k]);
    }
}

inline
void
CoordSys::HiFace (const IntVect& point,
                  int            dir,
                  Array<Real>&   loc) const
{
    loc.resize(BL_SPACEDIM);
    HiFace(point,dir, loc.dataPtr());
}

inline
void
CoordSys::LoNode (const IntVect& point,
                  Real*          loc) const
{
    assert(ok);
    assert(loc != 0);
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        loc[k] = offset[k] + dx[k]*point[k];
    }
}

inline
void
CoordSys::LoNode (const IntVect& point,
                  Array<Real>&   loc) const
{
    loc.resize(BL_SPACEDIM);
    LoNode(point, loc.dataPtr());
}

inline
void
CoordSys::HiNode (const IntVect& point,
                  Real*          loc) const
{
    assert(ok);
    assert(loc != 0);
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        loc[k] = offset[k] + dx[k]*(point[k] + 1);
    }
}

inline
void
CoordSys::HiNode (const IntVect& point,
                  Array<Real>&   loc) const
{
    loc.resize(BL_SPACEDIM);
    HiNode(point, loc.dataPtr());
}

inline
IntVect
CoordSys::CellIndex (const Real* point) const
{
    assert(ok);
    assert(point != 0);
    IntVect ix;
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        ix[k] = (int) ((point[k]-offset[k])/dx[k]);
    }
    return ix;
}

inline
IntVect
CoordSys::LowerIndex (const Real* point) const
{
    assert(ok);
    assert(point != 0);
    IntVect ix;
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        ix[k] = (int) ((point[k]-offset[k])/dx[k]);
    }
    return ix;
}

inline
IntVect
CoordSys::UpperIndex(const Real* point) const
{
    assert(ok);
    assert(point != 0);
    IntVect ix;
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        ix[k] = (int) ((point[k]-offset[k])/dx[k]);
    }
    return ix;
}

#endif /*_COORDSYS_H_*/
