//BL_COPYRIGHT_NOTICE

#ifndef _COORDSYS_H_
#define _COORDSYS_H_ 

//
// $Id: CoordSys.H,v 1.1 1997-12-10 19:07:40 lijewski Exp $
//

#include <REAL.H>
#include <Array.H>
#include <Box.H>
#include <FArrayBox.H>

class CoordSys
{
public:
#if (BL_SPACEDIM == 2)    
    enum CoordType { undef = -1, cartesian = 0, RZ = 1 };
#endif
#if (BL_SPACEDIM == 3)
    enum CoordType { undef = -1, cartesian = 0 };
#endif    

      // -------------------------------------------
      // ---------  I/O friend functions  ----------
      // -------------------------------------------
    friend ostream& operator << (ostream&, const CoordSys& );
    friend istream& operator >> (istream&, CoordSys& );

      // Default constructor, CoordSys is undefined
    CoordSys();

      // specify a cell size, and low end of index range
    CoordSys(const Real* cell_dx);
    
    ~CoordSys() {};

      // initialize after construction
    void define(const Real* cell_dx);

      // functions to set and access static members
    static void SetCoord( CoordType coord );
    static CoordType Coord() {
	return c_sys;
    }

    static int IsRZ();
    static int IsCartesian();

    static void SetOffset( const Real* x_lo );
    static const Real* Offset() {
	return offset;
    }
    static Real Offset(int dir) {
	return offset[dir];
    }

      // functions to access non-static members
    const Real* CellSize() const {
	assert(ok);
	return dx;
    }
    Real CellSize(int dir) const {
	assert(ok);
	return dx[dir];
    }

      // Index space -> Physical space mappings
    void CellCenter(const IntVect& point, Array<Real>& loc) const;
    void CellCenter(const IntVect& point, Real *loc) const;
    Real CellCenter(int point, int dir) const;
    Real LoEdge(int point, int dir) const;
    Real LoEdge(const IntVect& point, int dir) const;
    Real HiEdge(int point, int dir) const;
    Real HiEdge(const IntVect& point, int dir) const;
      // location of lo and hi cell faces
    void LoFace(const IntVect& point, int dir, Array<Real>& loc) const;
    void LoFace(const IntVect& point, int dir, Real *loc) const;
    void HiFace(const IntVect& point, int dir, Array<Real>& loc) const;
    void HiFace(const IntVect& point, int dir, Real *loc) const;
      // return location of lower left hand corner
    void LoNode(const IntVect& point, Array<Real>& loc) const;
    void LoNode(const IntVect& point, Real *loc) const;
      // return location of upper right hand corner
    void HiNode(const IntVect& point, Array<Real>& loc) const;
    void HiNode(const IntVect& point, Real *loc) const;

      // Physical Space -> index space mappings
      // WARNING:: May return undesired results if point
      //           is on a cell boundary.
      // returns cell centered index of cell containing point
    IntVect CellIndex(const Real* point) const;
      // returns node centered index of lower left hand corner of
      // cell containing this point
    IntVect LowerIndex(const Real* point) const;
      // returns node centered index of upper right hand corner of
      // cell containing this point
    IntVect UpperIndex(const Real* point) const;
    
      // compute cell volumes in given region
      // this version resizes input FAB
    void GetVolume( FArrayBox& vol, const Box& region) const;

      // this version returns volume in new FAB
    FArrayBox* GetVolume( const Box& region) const;
    
#if (BL_SPACEDIM == 2)
  //@ManDoc: compute d(log(A))/dr at cell centers in given region -- this version resizes input FAB
  void GetDLogA( FArrayBox& dloga, const Box& region, int dir) const;

  //@ManDoc: this version returns dloga in new FAB
  FArrayBox* GetDLogA( const Box& region, int dir) const;
#endif

      // these versions return volume of single cell
    Real Volume(const IntVect& point) const;
    Real Volume(const Real xlo[BL_SPACEDIM], const Real xhi[BL_SPACEDIM]) const;

      // compute area of cell faces in given region and given
      // index direction
    void GetFaceArea( FArrayBox& area, const Box& region, int dir) const;
    FArrayBox* GetFaceArea(const Box& region, int dir) const;
    
      // these versions return face areas of given cell in direction dir
    Real AreaLo(const IntVect& point, int dir) const;
    Real AreaHi(const IntVect& point, int dir) const;
    
      // return array of physical locations of cell edges
      // resizes input array
    void GetEdgeLoc( Array<Real>& loc,
		     const Box& region, int dir) const;

      // return array of physical locations of cell centers
      // resizes input array
    void GetCellLoc( Array<Real>& loc,
		     const Box& region, int dir) const;

      // return array of volume coordinates at cell edges
      // for region in given direction
    void GetEdgeVolCoord( Array<Real>& loc,
			  const Box& region, int dir) const;

      // return array of volume coordinates at cell centers
      // for region in given direction
    void GetCellVolCoord( Array<Real>& loc,
			  const Box& region, int dir) const;

protected:
    static CoordType c_sys;
    static Real  offset[BL_SPACEDIM];
        
    int        ok;
    Real       dx[BL_SPACEDIM];
};

#endif

