//BL_COPYRIGHT_NOTICE

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

//
// $Id: Geometry.H,v 1.1 1997-12-10 19:07:41 lijewski Exp $
//

#include <CoordSys.H>
#include <RealBox.H>
#include <BoxArray.H>
#include <MultiFab.H>

// This class describes problem domain and coordinate system for
// RECTANGULAR problem domains.  Since the problem domain is RECTANGULAR,
// periodicity is meaningful

class Geometry : public CoordSys
{
protected:
  static bool is_periodic[BL_SPACEDIM];// 0 means not periodic
    static RealBox prob_domain;
    Box    domain;

    void   read_params();
    
public:
#ifdef __GNUG__
     bool operator== (const Geometry& rhs) const { return true; }
#endif 
      // -------------------------------------------
      // ---------  I/O friend functions  ----------
      // -------------------------------------------
    friend ostream& operator << (ostream&, const Geometry& );
    friend istream& operator >> (istream&, Geometry& );

    Geometry();
    Geometry(const Box& dom);
    Geometry(const Geometry& g);
    
    ~Geometry() {};

      // read static values from ParmParse database.
    static void Setup();

    void define(const Box& dom);

    static const RealBox& ProbDomain()
    {
	return prob_domain;
    }
    
    static const Real* ProbLo();
    static const Real* ProbHi();
    static Real ProbLo(int dir);
    static Real ProbHi(int dir);
    static const Real* ProbLength();
    static Real ProbLength(int dir);

    const Box& Domain() const {
	return domain;
    }

      // define a multifab of areas and volumes with given grow factor
    void GetVolume(MultiFab& vol, const BoxArray& grds, int grow) const;
#if (BL_SPACEDIM == 2)
    void GetDLogA(MultiFab& dloga, const BoxArray& grds, int dir, 
		  int grow) const;
#endif
    void GetFaceArea(MultiFab& area, const BoxArray& grds, int dir,
		     int grow) const;

  // periodic info
  static bool isPeriodic( int dir ){
    return is_periodic[dir]!=0;
  }
  static bool isAnyPeriodic( ){
    return isPeriodic(0)
#if BL_SPACEDIM>1
      ||   isPeriodic(1)
#endif
#if BL_SPACEDIM>2
      ||   isPeriodic(2)
#endif
      ;
  }
  int period( int dir ) {
    assert( is_periodic[dir] );
    return domain.length(dir);
  }
  // compute Array of shifts which will translate src so that it will
  // intersect target with non-zero intersection.  the array will be
  // resized internally, so anything previously there will be gone
  // DO NOT return non-periodic shifts, even if the box's do intersect
  // without shifting.  The logic is that you will only do this as
  // a special case if there is some periodicity.  
  void periodicShift( const Box& target, const Box& src, 
		      Array<IntVect>& out ) const;
	
};

#endif
