//BL_COPYRIGHT_NOTICE

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

//
// $Id: Geometry.H,v 1.31 1999-04-13 15:28:19 lijewski Exp $
//

#include <CoordSys.H>
#include <RealBox.H>
#include <MultiFab.H>

#ifdef BL_USE_NEW_HFILES
#include <list>
using std::ostream;
using std::list;
#else
#include <list.h>
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
    explicit Geometry (const Box& dom);
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
    //@ManDoc: Fill ghost cells of all components with periodic data.
    //
    void FillPeriodicBoundary (MultiFab& mf,
                               bool      no_overlap,
                               bool      do_corners) const;
    //
    //@ManDoc: Fill ghost cells of selected components with periodic data.
    //
    void FillPeriodicBoundary (MultiFab& mf,
                               int       src_comp,
                               int       num_comp,
                               bool      no_overlap,
                               bool      do_corners) const;
    //
    // Src/dest box pairs to copy for periodic.
    //
    struct PIRec
    {
	PIRec ()
            :
            mfid(-1),
            srcId(-1) {}

	PIRec (int        _mfid,
               int        _srcId,
               const Box& _srcBox,
               const Box& _dstBox)
	    :
            mfid(_mfid),
            srcId(_srcId),
            srcBox(_srcBox),
            dstBox(_dstBox) {}

        int       mfid;
	int       srcId;
	Box       srcBox;
	Box       dstBox;
	FillBoxId fbid;

	friend ostream& operator<< (ostream&      os,
                                    const PIRec&  pir);
    };
    //
    //@ManDoc: Handy typedef.
    //
    typedef vector<PIRec> PIRMMap;
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
    //@ManDoc: The size of the PIRM cache.
    //
    static int PIRMCacheSize ();

protected:
    //
    // Used in caching PIRMMaps && CommData.
    //
    struct FPB
    {
        FPB ();

        FPB (const BoxArray& ba,
             const Box&      domain,
             int             scomp,
             int             ncomp,
             int             ngrow,
             bool            no_overlap,
             bool            do_corners);

        FPB (const FPB& rhs);

        ~FPB ();

        bool operator== (const FPB& rhs) const;
        bool operator!= (const FPB& rhs) const;

        Array<int>    m_cache;    // Snds cached for CollectData().
        CommDataCache m_commdata; // Yet another cache for CollectData().
        PIRMMap       m_pirm;
        BoxArray      m_ba;
        Box           m_domain;
        int           m_scomp;
        int           m_ncomp;
        int           m_ngrow;
        bool          m_no_overlap;
        bool          m_do_corners;
    };
    //
    // A useful typedef.
    //
    typedef list<FPB> FPBList;
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
    static bool    is_periodic[BL_SPACEDIM]; // 0 means not periodic
    static RealBox prob_domain;
    static FPBList m_FPBCache;
    //
    // Non-static data.
    //
    Box domain;
};

//
// Inlines.
//

inline
Geometry::FPB::FPB ()
    :
    m_scomp(-1),
    m_ncomp(-1),
    m_ngrow(-1),
    m_no_overlap(false),
    m_do_corners(false)
{}

inline
Geometry::FPB::FPB (const BoxArray& ba,
                    const Box&      domain,
                    int             scomp,
                    int             ncomp,
                    int             ngrow,
                    bool            no_overlap,
                    bool            do_corners)
    :
    m_ba(ba),
    m_domain(domain),
    m_scomp(scomp),
    m_ncomp(ncomp),
    m_ngrow(ngrow),
    m_no_overlap(no_overlap),
    m_do_corners(do_corners)
{
    assert(scomp >= 0);
    assert(ncomp >  0);
    assert(ngrow >= 0);
    assert(domain.ok());
}

inline
Geometry::FPB::FPB (const FPB& rhs)
    :
    m_cache(rhs.m_cache),
    m_commdata(rhs.m_commdata),
    m_pirm(rhs.m_pirm),
    m_ba(rhs.m_ba),
    m_domain(rhs.m_domain),
    m_scomp(rhs.m_scomp),
    m_ncomp(rhs.m_ncomp),
    m_ngrow(rhs.m_ngrow),
    m_no_overlap(rhs.m_no_overlap),
    m_do_corners(rhs.m_do_corners)
{}

inline
Geometry::FPB::~FPB () {}

inline
bool
Geometry::FPB::operator== (const FPB& rhs) const
{
    return
        m_scomp      == rhs.m_scomp      &&
        m_ncomp      == rhs.m_ncomp      &&
        m_ngrow      == rhs.m_ngrow      &&
        m_no_overlap == rhs.m_no_overlap &&
        m_do_corners == rhs.m_do_corners &&
        m_domain     == rhs.m_domain     &&
        m_ba         == rhs.m_ba;
}

inline
bool
Geometry::FPB::operator!= (const FPB& rhs) const
{
    return !operator==(rhs);
}

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
Geometry::FillPeriodicBoundary (MultiFab& mf,
                                bool no_overlap,
                                bool do_corners) const
{
    FillPeriodicBoundary(mf, 0, mf.nComp(), no_overlap, do_corners);
}

inline
Geometry::FPB&
Geometry::getFPB (MultiFab&  mf,
                  const FPB& fpb) const
{
    assert(isAnyPeriodic());
    //
    // Have we already made one with appropriate characteristics?
    //
    // Search from the front to the back ...
    //
    FPBList::iterator it = m_FPBCache.begin();

    for ( ; it != m_FPBCache.end(); ++it)
        if (*it == fpb)
            return *it;

    return buildFPB(mf,fpb);
}

#endif /*_GEOMETRY_H_*/
