//BL_COPYRIGHT_NOTICE

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

//
// $Id: Geometry.H,v 1.2 1997-12-11 05:01:09 lijewski Exp $
//

#include <CoordSys.H>
#include <RealBox.H>

//
// Forward declarations.
//
class BoxArray;
class MultiFab;

//
// This class describes problem domain and coordinate system for
// RECTANGULAR problem domains.  Since the problem domain is RECTANGULAR,
// periodicity is meaningful.
//

class Geometry
    :
    public CoordSys
{
public:

#ifdef __GNUG__
     bool operator== (const Geometry& rhs) const { return true; }
#endif 

    friend ostream& operator<< (ostream&, const Geometry& );
    friend istream& operator>> (istream&, Geometry& );

    Geometry ();
    Geometry (const Box& dom);
    Geometry (const Geometry& g);
    
    ~Geometry();
    //
    // Read static values from ParmParse database.
    //
    static void Setup ();

    void define (const Box& dom);

    static const RealBox& ProbDomain ();

    static const Real* ProbLo ();
    static const Real* ProbHi ();
    static Real ProbLo (int dir);
    static Real ProbHi (int dir);
    static const Real* ProbLength ();
    static Real ProbLength (int dir);

    const Box& Domain () const;
    //
    // Define a multifab of areas and volumes with given grow factor.
    //
    void GetVolume (MultiFab&       vol,
                    const BoxArray& grds,
                    int             grow) const;

#if (BL_SPACEDIM == 2)
    void GetDLogA (MultiFab&       dloga,
                   const BoxArray& grds,
                   int             dir,
                   int             grow) const;
#endif
    void GetFaceArea (MultiFab&       area,
                      const BoxArray& grds,
                      int             dir,
                      int             grow) const;
    //
    // Periodic info.
    //
    static bool isPeriodic (int dir);

    static bool isAnyPeriodic ();

    int period (int dir);
    //
    // compute Array of shifts which will translate src so that it will
    // intersect target with non-zero intersection.  the array will be
    // resized internally, so anything previously there will be gone
    // DO NOT return non-periodic shifts, even if the box's do intersect
    // without shifting.  The logic is that you will only do this as
    // a special case if there is some periodicity.
    //
    void periodicShift (const Box&      target,
                        const Box&      src, 
                        Array<IntVect>& out) const;
protected:
    //
    // Helper function.
    //
    void read_params ();
    //
    // Static data.
    //
    static bool is_periodic[BL_SPACEDIM]; // 0 means not periodic
    static RealBox prob_domain;
    //
    // Non-static data.
    //
    Box domain;
};

//
// Inlines.
//
inline
Geometry::~Geometry() {}

inline
const RealBox&
Geometry::ProbDomain ()
{
    return prob_domain;
}

inline
const Box&
Geometry::Domain () const
{
    return domain;
}

inline
bool
Geometry::isPeriodic (int dir)
{
    return is_periodic[dir]!=0;
}

inline
bool
Geometry::isAnyPeriodic ()
{
    return isPeriodic(0)
#if BL_SPACEDIM>1
        ||   isPeriodic(1)
#endif
#if BL_SPACEDIM>2
        ||   isPeriodic(2)
#endif
        ;
}

inline
int
Geometry::period (int dir)
{
    assert(is_periodic[dir]);
    return domain.length(dir);
}

inline
const Real*
Geometry::ProbLo()
{
    return prob_domain.lo();
}

inline
const Real*
Geometry::ProbHi()
{
    return prob_domain.hi();
}

inline
Real
Geometry::ProbLo(int dir)
{
    return prob_domain.lo(dir);
}

inline
Real
Geometry::ProbHi(int dir)
{
    return prob_domain.hi(dir);
}

inline
const Real*
Geometry::ProbLength()
{
    return prob_domain.length();
}

inline
Real
Geometry::ProbLength(int dir)
{
    return prob_domain.length(dir);
}

inline
Geometry::Geometry () {}

inline
Geometry::Geometry (const Box& dom)
{
    define(dom);
}

inline
Geometry::Geometry (const Geometry& g)
{
    define(g.domain);
}

#endif /*_GEOMETRY_H_*/
