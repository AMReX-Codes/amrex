//BL_COPYRIGHT_NOTICE

#ifndef _COORDSYS_H_
#define _COORDSYS_H_ 

//
// $Id: CoordSys.H,v 1.3 1997-12-11 23:26:59 lijewski Exp $
//

#include <REAL.H>
#include <Array.H>
#include <Box.H>

//
// Forward declaration.
//
class FArrayBox;

class CoordSys
{
public:
#if (BL_SPACEDIM == 2)    
    enum CoordType { undef = -1, cartesian = 0, RZ = 1 };
#endif
#if (BL_SPACEDIM == 3)
    enum CoordType { undef = -1, cartesian = 0 };
#endif    

    friend ostream& operator<< (ostream&, const CoordSys& );
    friend istream& operator>> (istream&, CoordSys& );
    //
    // Default constructor, CoordSys is undefined.
    //
    CoordSys ();
    //
    // Specify a cell size, and low end of index range.
    //
    CoordSys (const Real* cell_dx);
    
    ~CoordSys ();
    //
    // Initialize after construction.
    //
    void define (const Real* cell_dx);

    static void SetCoord (CoordType coord);

    static CoordType Coord ();

    static bool IsRZ ();
    static bool IsCartesian ();

    static void SetOffset (const Real* x_lo);
    static const Real* Offset ();

    static Real Offset(int dir);
    //
    // Functions to access non-static members.
    //
    const Real* CellSize() const;

    Real CellSize(int dir) const;
    //
    // Index space -> Physical space mappings.
    //
    void CellCenter (const IntVect& point, Array<Real>& loc) const;
    void CellCenter (const IntVect& point, Real *loc) const;
    Real CellCenter (int point, int dir) const;
    Real LoEdge (int point, int dir) const;
    Real LoEdge (const IntVect& point, int dir) const;
    Real HiEdge (int point, int dir) const;
    Real HiEdge (const IntVect& point, int dir) const;
    //
    // Location of lo and hi cell faces.
    //
    void LoFace (const IntVect& point, int dir, Array<Real>& loc) const;
    void LoFace (const IntVect& point, int dir, Real *loc) const;
    void HiFace (const IntVect& point, int dir, Array<Real>& loc) const;
    void HiFace (const IntVect& point, int dir, Real *loc) const;
    //
    // Return location of lower left hand corner.
    //
    void LoNode (const IntVect& point, Array<Real>& loc) const;
    void LoNode (const IntVect& point, Real *loc) const;
    //
    // Return location of upper right hand corner.
    //
    void HiNode (const IntVect& point, Array<Real>& loc) const;
    void HiNode (const IntVect& point, Real *loc) const;
    //
    // Physical Space -> index space mappings
    // WARNING:: May return undesired results if point
    //           is on a cell boundary.
    // returns cell centered index of cell containing point.
    //
    IntVect CellIndex (const Real* point) const;
    //
    // Returns node centered index of lower left hand corner of
    // cell containing this point.
    //
    IntVect LowerIndex (const Real* point) const;
    //
    // returns node centered index of upper right hand corner of
    // cell containing this point.
    //
    IntVect UpperIndex (const Real* point) const;
    //
    // Compute cell volumes in given region this version resizes input FAB.
    //
    void GetVolume (FArrayBox& vol, const Box& region) const;
    //
    // This version returns volume in new FAB.
    //
    FArrayBox* GetVolume (const Box& region) const;
    
#if (BL_SPACEDIM == 2)
    //
    // Compute d(log(A))/dr at cell centers in given region
    // This version resizes input FAB.
    //
    void GetDLogA (FArrayBox& dloga, const Box& region, int dir) const;
    //
    // This version returns dloga in new FAB.
    //
    FArrayBox* GetDLogA( const Box& region, int dir) const;
#endif
    //
    // These versions return volume of single cell.
    //
    Real Volume (const IntVect& point) const;
    Real Volume (const Real xlo[BL_SPACEDIM],
                 const Real xhi[BL_SPACEDIM]) const;
    //
    // compute area of cell faces in given region and given index direction.
    //
    void GetFaceArea (FArrayBox& area, const Box& region, int dir) const;
    FArrayBox* GetFaceArea (const Box& region, int dir) const;
    //
    // These versions return face areas of given cell in direction dir.
    //
    Real AreaLo (const IntVect& point, int dir) const;
    Real AreaHi (const IntVect& point, int dir) const;
    //
    // Return array of physical locations of cell edges resizes input array.
    //
    void GetEdgeLoc (Array<Real>& loc, const Box& region, int dir) const;
    //
    // Return array of physical locations of cell centers resizes input array.
    //
    void GetCellLoc (Array<Real>& loc, const Box& region, int dir) const;
    //
    // Return array of volume coordinates at cell edges
    // for region in given direction.
    //
    void GetEdgeVolCoord (Array<Real>& loc, const Box& region, int dir) const;
    //
    // Return array of volume coordinates at cell centers
    // for region in given direction.
    //
    void GetCellVolCoord (Array<Real>& loc, const Box& region, int dir) const;

protected:
    //
    // Static data.
    //
    static CoordType c_sys;
    static Real  offset[BL_SPACEDIM];
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
    CellCenter(point,&(loc[0]));
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
    LoFace(point,dir,&(loc[0]));
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
    HiFace(point,dir,&(loc[0]));
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
    LoNode(point,&(loc[0]));
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
    HiNode(point,&(loc[0]));
}

#endif /*_COORDSYS_H_*/
