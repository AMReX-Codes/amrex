#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_
//
// $Id: Geometry.H,v 1.40 2002-04-27 07:07:41 marc Exp $
//
#include <iosfwd>
#include <list>

#include <CoordSys.H>
#include <MultiFab.H>
#include <ParallelDescriptor.H>
#include <RealBox.H>

struct FPB;

//
//@Man:
//@Memo: Rectangular problem domain geometry.
/*@Doc:

  This class describes problem domain and coordinate system for
  RECTANGULAR problem domains.  Since the problem domain is RECTANGULAR,
  periodicity is meaningful.
*/

class Geometry
    :
    public CoordSys
{
public:
    //
    //@ManDoc: The default constructor.
    //
    Geometry ();
    //
    //@ManDoc: Constructor taking the rectangular domain.
    //
    //explicit Geometry (const Box& dom);
    Geometry (const Box& dom, const RealBox* rb = 0, int coord = -1);
    //
    //@ManDoc: The copy constructor.
    //
    Geometry (const Geometry& g);
    //
    //@ManDoc: The destructor.
    //
    ~Geometry();
    //
    //@ManDoc: Read static values from ParmParse database.
    //
    static void Setup (const RealBox* rb = 0, int coord = -1);
    //
    //@ManDoc: Set the rectangular domain after using default constructor.
    //
    void define (const Box& dom, const RealBox* rb = 0, int coord = -1);
    //
    //@ManDoc: Returns the problem domain.
    //
    static const RealBox& ProbDomain ();
    //
    //@ManDoc: Sets the problem domain.
    //
    static void ProbDomain (const RealBox& rb);
    //
    //@ManDoc: Returns the lo end of the problem domain in each dimension.
    //
    static const Real* ProbLo ();
    //
    //@ManDoc: Returns the hi end of the problem domain in each dimension.
    //
    static const Real* ProbHi ();
    //
    //@ManDoc: Returns the lo end of the problem domain in specified direction.
    //
    static Real ProbLo (int dir);
    //
    //@ManDoc: Returns the hi end of the problem domain in specified direction.
    //
    static Real ProbHi (int dir);
    //
    //@ManDoc: Returns length of problem domain in each dimension.
    //
    static const Real* ProbLength ();
    //
    //@ManDoc: Returns length of problem domain in specified dimension.
    //
    static Real ProbLength (int dir);
    //
    //@ManDoc: Returns our rectangular domain.
    //
    const Box& Domain () const;
    //
    //@ManDoc: Sets our rectangular domain.
    //
    void Domain (const Box& bx);
    //
    //@ManDoc: Define a multifab of areas and volumes with given grow factor.
    //
    void GetVolume (MultiFab&       vol,
                    const BoxArray& grds,
                    int             grow) const;

    /*@ManDoc: Compute d(log(A))/dr at cell centers in given region and
               stuff the results into the passed MultiFab.
    */
    void GetDLogA (MultiFab&       dloga,
                   const BoxArray& grds,
                   int             dir,
                   int             grow) const;

    /*@ManDoc: Compute area of cell faces in given region and stuff
               stuff the results into the passed MultiFab.
    */
    void GetFaceArea (MultiFab&       area,
                      const BoxArray& grds,
                      int             dir,
                      int             grow) const;
    //
    //@ManDoc: Is the domain periodic in the specified direction?
    //
    static bool isPeriodic (int dir);
    //
    //@ManDoc: Is domain periodic in any direction?
    //
    static bool isAnyPeriodic ();
    //
    //@ManDoc: What's period in specified direction?
    //
    int period (int dir) const;

    /*@ManDoc: Compute Array of shifts which will translate src so that it will
               intersect target with non-zero intersection.  the array will be
               resized internally, so anything previously there will be gone
               DO NOT return non-periodic shifts, even if the box's do
               intersect without shifting.  The logic is that you will only do
               this as a special case if there is some periodicity.
    */
    void periodicShift (const Box&      target,
                        const Box&      src, 
                        Array<IntVect>& out) const;
    //
    //@ManDoc: Fill ghost cells of all components with periodic data.
    //
    void FillPeriodicBoundary (MultiFab& mf,
                               bool      do_corners = false) const;
    //
    //@ManDoc: Fill ghost cells of selected components with periodic data.
    //
    void FillPeriodicBoundary (MultiFab& mf,
                               int       src_comp,
                               int       num_comp,
                               bool      do_corners = false) const;
    //
    //@ManDoc: Flush the cache of PIRM information.
    //
    static void FlushPIRMCache ();
    //
    //@ManDoc: The size of the PIRM cache.
    //
    static int PIRMCacheSize ();

protected:
    //
    // Helper functions.
    //
    void read_params ();
    //
    // Helper functions for FillPeriodicBoundary(). 
    //
    FPB& getFPB (MultiFab&  mf,
                 const FPB& fpb) const;

    FPB& buildFPB (MultiFab&  mf,
                   const FPB& fpb) const;
    //
    // Static data.
    //
    static int     spherical_origin_fix;
    static bool    is_periodic[BL_SPACEDIM]; // 0 means not periodic
    static RealBox prob_domain;
    //
    // Non-static data.
    //
    Box domain;
};
//
//@ManDoc: Nice ASCII output.
//
std::ostream& operator<< (std::ostream&, const Geometry&);
//
//@ManDoc: Nice ASCII input.
//
std::istream& operator>> (std::istream&, Geometry&);

#endif /*_GEOMETRY_H_*/
