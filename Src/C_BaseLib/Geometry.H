//BL_COPYRIGHT_NOTICE

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

//
// $Id: Geometry.H,v 1.14 1998-06-13 17:28:30 lijewski Exp $
//

#include <CoordSys.H>
#include <RealBox.H>
#include <MultiFab.H>

#ifdef BL_USE_NEW_HFILES
#include <list>
#include <map>
using std::ostream;
using std::multimap;
using std::list;
#else
#include <list.h>
#include <map.h>
#endif

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

#ifdef __GNUG__
    //
    // Hack to work around g++ bug :-(
    //
    bool operator== (const Geometry& rhs) const { return true; }
#endif 
    //
    //@ManDoc: Nice ASCII output.
    //
    friend ostream& operator<< (ostream&, const Geometry& );
    //
    //@ManDoc: Nice ASCII input.
    //
    friend istream& operator>> (istream&, Geometry& );
    //
    //@ManDoc: The default constructor.
    //
    Geometry ();
    //
    //@ManDoc: Constructor taking the rectangular domain.
    //
    Geometry (const Box& dom);
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
    static void Setup ();
    //
    //@ManDoc: Set the rectangular domain after using default constructor.
    //
    void define (const Box& dom);
    //
    //@ManDoc: Returns the problem domain.
    //
    static const RealBox& ProbDomain ();
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
    //@ManDoc: Define a multifab of areas and volumes with given grow factor.
    //
    void GetVolume (MultiFab&       vol,
                    const BoxArray& grds,
                    int             grow) const;

#if (BL_SPACEDIM == 2)
    /*@ManDoc: Compute d(log(A))/dr at cell centers in given region and
               stuff the results into the passed MultiFab.
    */
    void GetDLogA (MultiFab&       dloga,
                   const BoxArray& grds,
                   int             dir,
                   int             grow) const;
#endif

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
    int period (int dir);

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
    //@ManDoc: Fill ghost cells with periodic data.
    //
    void FillPeriodicBoundary (MultiFab& mf,
                               bool      no_ovlp = false) const;
    //
    //@ManDoc: Fill ghost cells with periodic data.
    //
    void FillPeriodicBoundary (MultiFab& mf,
                               int       src_comp,
                               int       num_comp,
                               bool      no_ovlp = false) const;
    //
    // Src/dest box pairs to copy for periodic.
    //
    class PIRec
    {
    public:
	PIRec ()
            :
            srcId(-1) {}
	PIRec (int        _srcId,
               const Box& _srcBox,
               const Box& _dstBox)
	    :
            srcId(_srcId),
            srcBox(_srcBox),
            dstBox(_dstBox) {}

	int       srcId;
	Box       srcBox;
	Box       dstBox;
	FillBoxId fbid;
#ifdef __GNUG__
	bool operator== (const PIRec& rhs) const { return true; }
#endif
	friend ostream& operator<< (ostream&      os,
                                    const PIRec&  pir);
    };
    //
    //@ManDoc: Handy typedef.
    //
    typedef multimap<int,PIRec> PIRMMap;
    //
    //@ManDoc: Handy output function.
    //
    friend ostream& operator<< (ostream&       os,
                                const PIRMMap& pirm);
    //
    //@ManDoc: Flush the cache of PIRM information.
    //
    static void FlushPIRMCache ();
    //
    //@ManDoc: The size of the cache.
    //
    static int PIRMCacheSize ();

    /*@ManDoc: Compute periodic intersection records for data into each
               fab on this proc.  It's an error to call this if
               !isAnyPeriodic().
    */
    PIRMMap& getPIRMMap (const BoxArray& grids,
                         int             nGrow,
                         bool            no_ovlp) const;
protected:
    //
    // Used in caching FillPeriodicBoundary information.
    //
    struct FPB
    {
        FPB ()
            :
            m_pirm(0),
            m_ngrow(-1),
            m_noovlp(false) {}

        FPB (const BoxArray& ba,
             int             ngrow,
             bool            noovlp)
            :
            m_pirm(0),
            m_ba(ba),
            m_ngrow(ngrow),
            m_noovlp(true)
        {
            assert(ngrow >= 0);
        }

        FPB (PIRMMap*        pirm,
             const BoxArray& ba,
             int             ngrow,
             bool            noovlp)
            :
            m_pirm(pirm),
            m_ba(ba),
            m_ngrow(ngrow),
            m_noovlp(true)
        {
            assert(ngrow >= 0);
            assert(!(pirm == 0));
        }

        ~FPB () { /* Don't delete m_pirm in destructor! */ }

        bool operator== (const FPB& rhs) const;

        PIRMMap* m_pirm;
        BoxArray m_ba;
        int      m_ngrow;
        bool     m_noovlp;
    };
    //
    // Helper functions.
    //
    void read_params ();

    //
    // Helper function for getPIRMMap().
    //
    PIRMMap& buildPIRMMap (const BoxArray& grids,
                           int             nGrow,
                           bool            no_ovlp) const;
    //
    // Static data.
    //
    static bool      is_periodic[BL_SPACEDIM]; // 0 means not periodic
    static RealBox   prob_domain;
    static list<FPB> m_Cache;
    //
    // Non-static data.
    //
    Box domain;
};

//
// Inlines.
//

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
Geometry::ProbLo ()
{
    return prob_domain.lo();
}

inline
const Real*
Geometry::ProbHi ()
{
    return prob_domain.hi();
}

inline
Real
Geometry::ProbLo (int dir)
{
    return prob_domain.lo(dir);
}

inline
Real
Geometry::ProbHi (int dir)
{
    return prob_domain.hi(dir);
}

inline
const Real*
Geometry::ProbLength ()
{
    return prob_domain.length();
}

inline
Real
Geometry::ProbLength (int dir)
{
    return prob_domain.length(dir);
}

inline
void
Geometry::FillPeriodicBoundary (MultiFab& mf, bool no_ovlp) const
{
    FillPeriodicBoundary(mf, 0, mf.nComp(), no_ovlp);
}

inline
int
Geometry::PIRMCacheSize ()
{
    return m_Cache.size();
}

#endif /*_GEOMETRY_H_*/
