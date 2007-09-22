//
// $Id: MultiFab.cpp,v 1.83 2007-09-22 15:05:31 lijewski Exp $
//
#include <winstd.H>

#include <algorithm>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <map>
#include <limits>

#include <BLassert.H>
#include <MultiFab.H>
#include <ParallelDescriptor.H>
#include <Profiler.H>
#include <ParmParse.H>

bool MultiFab::check_for_nan = false;
bool MultiFab::check_for_inf = false;

void
MultiFab::Add (MultiFab&       dst,
	       const MultiFab& src,
	       int             srccomp,
	       int             dstcomp,
	       int             numcomp,
	       int             nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrow() >= nghost && src.nGrow() >= nghost);

    for (MFIter mfi(dst); mfi.isValid(); ++mfi)
    {
        Box bx = BoxLib::grow(mfi.validbox(),nghost);

        if (bx.ok())
            dst[mfi].plus(src[mfi], bx, bx, srccomp, dstcomp, numcomp);
    }
}

void
MultiFab::Copy (MultiFab&       dst,
                const MultiFab& src,
                int             srccomp,
                int             dstcomp,
                int             numcomp,
                int             nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrow() >= nghost && src.nGrow() >= nghost);

    for (MFIter mfi(dst); mfi.isValid(); ++mfi)
    {
        Box bx = BoxLib::grow(mfi.validbox(),nghost);

        if (bx.ok())
            dst[mfi].copy(src[mfi], bx, srccomp, bx, dstcomp, numcomp);
    }
}

void
MultiFab::Subtract (MultiFab&       dst,
		    const MultiFab& src,
		    int             srccomp,
		    int             dstcomp,
		    int             numcomp,
		    int             nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrow() >= nghost && src.nGrow() >= nghost);

    for (MFIter mfi(dst); mfi.isValid(); ++mfi)
    {
        Box bx = BoxLib::grow(mfi.validbox(),nghost);

        if (bx.ok())
            dst[mfi].minus(src[mfi], bx, bx, srccomp, dstcomp, numcomp);
    }
}

void
MultiFab::FillBoundary (bool local)
{
    FillBoundary(0, n_comp, local);
}

void
MultiFab::plus (Real val,
                int  nghost)
{
    plus(val,0,n_comp,nghost);
}

void
MultiFab::plus (Real       val,
                const Box& region,
                int        nghost)
{
    plus(val,region,0,n_comp,nghost);
}

void
MultiFab::mult (Real val,
                int  nghost)
{
    mult(val,0,n_comp,nghost);
}

void
MultiFab::mult (Real       val,
                const Box& region,
                int        nghost)
{
    mult(val,region,0,n_comp,nghost);
}

void
MultiFab::invert (Real numerator,
                  int  nghost)
{
    invert(numerator,0,n_comp,nghost);
}

void
MultiFab::invert (Real       numerator,
                  const Box& region,
                  int        nghost)
{
    invert(numerator,region,0,n_comp,nghost);
}

void
MultiFab::negate (int nghost)
{
    negate(0,n_comp,nghost);
}

void
MultiFab::negate (const Box& region,
                  int        nghost)
{
    negate(region,0,n_comp,nghost);
}

void
MultiFab::startup ()
{
    static bool first = true;

    if (first)
    {
        first = false;
        ParmParse pp("multifab");
        pp.query("check_for_nan", check_for_nan);
        pp.query("check_for_inf", check_for_inf);
    }
}

MultiFab::MultiFab ()
{
    startup();
}

MultiFab::MultiFab (const BoxArray& bxs,
                    int             ncomp,
                    int             ngrow,
                    FabAlloc        alloc)
    :
    FabArray<FArrayBox>(bxs,ncomp,ngrow,alloc)
{
    startup();

    if ((check_for_nan || check_for_inf) && alloc == Fab_allocate) setVal(0);
}

void
MultiFab::operator= (const Real& r)
{
    setVal(r);
}

void
MultiFab::define (const BoxArray& bxs,
                  int             nvar,
                  int             ngrow,
                  FabAlloc        alloc)
{
    this->FabArray<FArrayBox>::define(bxs,nvar,ngrow,alloc);

    if ((check_for_nan || check_for_inf) && alloc == Fab_allocate) setVal(0);
}

void
MultiFab::define (const BoxArray&            bxs,
                  int                        nvar,
                  int                        ngrow,
                  const DistributionMapping& dm,
                  FabAlloc                   alloc)
{
    this->FabArray<FArrayBox>::define(bxs,nvar,ngrow,dm,alloc);

    if ((check_for_nan || check_for_inf) && alloc == Fab_allocate) setVal(0);
}

bool 
MultiFab::contains_nan () const
{
    bool r = false;

    for (MFIter mfi(*this); mfi.isValid() && !r; ++mfi)
        if (fabparray[mfi.index()].contains_nan())
            r = true;

    ParallelDescriptor::ReduceBoolOr(r);

    return r;
}

bool 
MultiFab::contains_inf () const
{
    bool r = false;

    for (MFIter mfi(*this); mfi.isValid() && !r; ++mfi)
        if (fabparray[mfi.index()].contains_inf())
            r = true;

    ParallelDescriptor::ReduceBoolOr(r);

    return r;
}

static
void
AbortOnNaN (const FArrayBox& fab)
{
    std::cout << fab << std::endl;
    BoxLib::Abort("FArrayBox contains a NaN");
}

static
void
AbortOnInf (const FArrayBox& fab)
{
    std::cout << fab << std::endl;
    BoxLib::Abort("FArrayBox contains a Inf");
}

const FArrayBox&
MultiFab::operator[] (const MFIter& mfi) const
{
    if (check_for_nan && fabparray[mfi.index()].contains_nan())
        AbortOnNaN(fabparray[mfi.index()]);

    if (check_for_inf && fabparray[mfi.index()].contains_inf())
        AbortOnInf(fabparray[mfi.index()]);

    return fabparray[mfi.index()];
}

const FArrayBox&
MultiFab::get (const MFIter& mfi) const
{
    if (check_for_nan && fabparray[mfi.index()].contains_nan())
        AbortOnNaN(fabparray[mfi.index()]);

    if (check_for_inf && fabparray[mfi.index()].contains_inf())
        AbortOnInf(fabparray[mfi.index()]);

    return fabparray[mfi.index()];
}

FArrayBox&
MultiFab::operator[] (const MFIter& mfi)
{
    if (check_for_nan && fabparray[mfi.index()].contains_nan())
        AbortOnNaN(fabparray[mfi.index()]);

    if (check_for_inf && fabparray[mfi.index()].contains_inf())
        AbortOnInf(fabparray[mfi.index()]);

    return fabparray[mfi.index()];
}

FArrayBox&
MultiFab::get (const MFIter& mfi)
{
    if (check_for_nan && fabparray[mfi.index()].contains_nan())
        AbortOnNaN(fabparray[mfi.index()]);

    if (check_for_inf && fabparray[mfi.index()].contains_inf())
        AbortOnInf(fabparray[mfi.index()]);

    return fabparray[mfi.index()];
}

const FArrayBox&
MultiFab::operator[] (int K) const
{
    if (check_for_nan && fabparray[K].contains_nan())
        AbortOnNaN(fabparray[K]);

    if (check_for_inf && fabparray[K].contains_inf())
        AbortOnInf(fabparray[K]);

    return fabparray[K];
}

FArrayBox&
MultiFab::operator[] (int K)
{
    if (check_for_nan && fabparray[K].contains_nan())
        AbortOnNaN(fabparray[K]);

    if (check_for_inf && fabparray[K].contains_inf())
        AbortOnInf(fabparray[K]);

    return fabparray[K];
}


Real
MultiFab::min (int comp,
               int nghost) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);

    Real mn = std::numeric_limits<Real>::max();

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        mn = std::min(mn,get(mfi).min(BoxLib::grow(mfi.validbox(),nghost),comp));
    }

    ParallelDescriptor::ReduceRealMin(mn);

    return mn;
}

Real
MultiFab::min (const Box& region,
               int        comp,
               int        nghost) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);

    Real mn = std::numeric_limits<Real>::max();

    for ( MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b = BoxLib::grow(mfi.validbox(),nghost) & region;

        if (b.ok())
            mn = std::min(mn, get(mfi).min(b,comp));
    }

    ParallelDescriptor::ReduceRealMin(mn);

    return mn;
}

Real
MultiFab::max (int comp,
               int nghost) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);

    Real mn = -std::numeric_limits<Real>::max();

    for ( MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        mn = std::max(mn, get(mfi).max(BoxLib::grow(mfi.validbox(),nghost),comp));
    }

    ParallelDescriptor::ReduceRealMax(mn);

    return mn;
}

Real
MultiFab::max (const Box& region,
               int        comp,
               int        nghost) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);

    Real mn = -std::numeric_limits<Real>::max();

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b = BoxLib::grow(mfi.validbox(),nghost) & region;

        if (b.ok())
            mn = std::max(mn, get(mfi).max(b,comp));
    }

    ParallelDescriptor::ReduceRealMax(mn);

    return mn;
}

void
MultiFab::minus (const MultiFab& mf,
                 int             strt_comp,
                 int             num_comp,
                 int             nghost)
{
    BL_ASSERT(boxarray == mf.boxarray);
    BL_ASSERT(strt_comp >= 0);
    BL_ASSERT(num_comp > 0);
#ifndef NDEBUG
    int lst_comp = strt_comp + num_comp - 1;
#endif
    BL_ASSERT(lst_comp < n_comp && lst_comp < mf.n_comp);
    BL_ASSERT(nghost <= n_grow && nghost <= mf.n_grow);

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        Box bx = BoxLib::grow(mfi.validbox(),nghost);

        get(mfi).minus(mf[mfi], bx, strt_comp, strt_comp, num_comp);
    }
}

void
MultiFab::divide (const MultiFab& mf,
		  int             strt_comp,
		  int             num_comp,
		  int             nghost)
{
    BL_ASSERT(boxarray == mf.boxarray);
    BL_ASSERT(strt_comp >= 0);
    BL_ASSERT(num_comp > 0);
#ifndef NDEBUG
    int lst_comp = strt_comp + num_comp - 1;
#endif
    BL_ASSERT(lst_comp < n_comp && lst_comp < mf.n_comp);
    BL_ASSERT(nghost <= n_grow && nghost <= mf.n_grow);

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        Box bx = BoxLib::grow(mfi.validbox(),nghost);

        get(mfi).divide(mf[mfi], bx, strt_comp, strt_comp, num_comp);
    }
}

void
MultiFab::plus (Real val,
                int  comp,
                int  num_comp,
                int  nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        get(mfi).plus(val,BoxLib::grow(mfi.validbox(),nghost),comp,num_comp);
    }
}

void
MultiFab::plus (Real       val,
                const Box& region,
                int        comp,
                int        num_comp,
                int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b = BoxLib::grow(mfi.validbox(),nghost) & region;

        if (b.ok())
            get(mfi).plus(val,b,comp,num_comp);
    }
}

void
MultiFab::plus (const MultiFab& mf,
                int             strt_comp,
                int             num_comp,
                int             nghost)
{
    BL_ASSERT(boxarray == mf.boxarray);
    BL_ASSERT(strt_comp >= 0);
    BL_ASSERT(num_comp > 0);
#ifndef NDEBUG
    int lst_comp = strt_comp + num_comp - 1;
#endif
    BL_ASSERT(lst_comp < n_comp && lst_comp < mf.n_comp);
    BL_ASSERT(nghost <= n_grow && nghost <= mf.n_grow);

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        Box bx = BoxLib::grow(mfi.validbox(),nghost);

        get(mfi).plus(mf[mfi], bx, strt_comp, strt_comp, num_comp);
    }
}

void
MultiFab::mult (Real val,
                int  comp,
                int  num_comp,
                int  nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        get(mfi).mult(val, BoxLib::grow(mfi.validbox(),nghost),comp,num_comp);
    }
}

void
MultiFab::mult (Real       val,
                const Box& region,
                int        comp,
                int        num_comp,
                int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b = BoxLib::grow(mfi.validbox(),nghost) & region;

        if (b.ok())
            get(mfi).mult(val, b, comp, num_comp);
    }
}

void
MultiFab::invert (Real numerator,
                  int  comp,
                  int  num_comp,
                  int  nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        get(mfi).invert(numerator, BoxLib::grow(mfi.validbox(),nghost),comp,num_comp);
    }
}

void
MultiFab::invert (Real       numerator,
                  const Box& region,
                  int        comp,
                  int        num_comp,
                  int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b = BoxLib::grow(mfi.validbox(),nghost) & region;

        if (b.ok())
            get(mfi).invert(numerator,b,comp,num_comp);
    }
}

void
MultiFab::negate (int comp,
                  int num_comp,
                  int nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+num_comp <= n_comp);

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        get(mfi).negate(BoxLib::grow(mfi.validbox(),nghost),comp,num_comp);
    }
}

void
MultiFab::negate (const Box& region,
                  int        comp,
                  int        num_comp,
                  int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+num_comp <= n_comp);

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b = BoxLib::grow(mfi.validbox(),nghost) & region;

        if (b.ok())
            get(mfi).negate(b,comp,num_comp);
    }
}

void
BoxLib::linInterpAddBox (MultiFabCopyDescriptor& fabCopyDesc,
                         BoxList*                returnUnfilledBoxes,
                         Array<FillBoxId>&       returnedFillBoxIds,
                         const Box&              subbox,
                         const MultiFabId&       faid1,
                         const MultiFabId&       faid2,
                         Real                    t1,
                         Real                    t2,
                         Real                    t,
                         int                     src_comp,
                         int                     dest_comp,
                         int                     num_comp,
                         bool                    extrap)
{
    const Real teps = (t2-t1)/1000.0;

    BL_ASSERT(extrap || ( (t>=t1-teps) && (t <= t2+teps) ) );

    if (t >= t1-teps && t <= t1+teps)
    {
        returnedFillBoxIds.resize(1);
        returnedFillBoxIds[0] = fabCopyDesc.AddBox(faid1,
                                                   subbox,
                                                   returnUnfilledBoxes,
                                                   src_comp,
                                                   dest_comp,
                                                   num_comp);
    }
    else if (t > t2-teps && t < t2+teps)
    {
        returnedFillBoxIds.resize(1);
        returnedFillBoxIds[0] = fabCopyDesc.AddBox(faid2,
                                                   subbox,
                                                   returnUnfilledBoxes,
                                                   src_comp,
                                                   dest_comp,
                                                   num_comp);
    }
    else
    {
        returnedFillBoxIds.resize(2);
        BoxList tempUnfilledBoxes(subbox.ixType());
        returnedFillBoxIds[0] = fabCopyDesc.AddBox(faid1,
                                                   subbox,
                                                   returnUnfilledBoxes,
                                                   src_comp,
                                                   dest_comp,
                                                   num_comp);
        returnedFillBoxIds[1] = fabCopyDesc.AddBox(faid2,
                                                   subbox,
                                                   &tempUnfilledBoxes,
                                                   src_comp,
                                                   dest_comp,
                                                   num_comp);
        //
        // The boxarrays for faid1 and faid2 should be the
        // same so only use returnUnfilledBoxes from one AddBox here.
        //
    }
}

void
BoxLib::linInterpFillFab (MultiFabCopyDescriptor& fabCopyDesc,
                          const Array<FillBoxId>& fillBoxIds,
                          const MultiFabId&       faid1,
                          const MultiFabId&       faid2,
                          FArrayBox&              dest,
                          Real                    t1,
                          Real                    t2,
                          Real                    t,
                          int                     src_comp,   // these comps need to be removed
                          int                     dest_comp,  // from this routine
                          int                     num_comp,
                          bool                    extrap)
{
    const Real teps = (t2-t1)/1000.0;

    BL_ASSERT(extrap || ( (t>=t1-teps) && (t <= t2+teps) ) );

    if (t >= t1-teps && t <= t1+teps)
    {
        fabCopyDesc.FillFab(faid1, fillBoxIds[0], dest);
    }
    else if (t > t2-teps && t < t2+teps)
    {
        fabCopyDesc.FillFab(faid2, fillBoxIds[0], dest);
    }
    else
    {
        BL_ASSERT(dest_comp + num_comp <= dest.nComp());

        FArrayBox dest1(dest.box(), dest.nComp());
        dest1.setVal(Real(1.e30)); // FIXME - Whats a better value?
        FArrayBox dest2(dest.box(), dest.nComp());
        dest2.setVal(Real(1.e30)); // FIXME - Whats a better value?
        fabCopyDesc.FillFab(faid1, fillBoxIds[0], dest1);
        fabCopyDesc.FillFab(faid2, fillBoxIds[1], dest2);
        dest.linInterp(dest1,
                       dest1.box(),
                       src_comp,
                       dest2,
                       dest2.box(),
                       src_comp,
                       t1,
                       t2,
                       t,
                       dest.box(),
                       dest_comp,
                       num_comp);
    }
}

//
// Holds single grid intersection record.
//

struct SIRec
{
    SIRec ()
        :
        m_i(-1),
        m_j(-1) {}

    SIRec (int        i,
           int        j,
           const Box& bx)
        :
        m_i(i),
        m_j(j),
        m_bx(bx)
    {
        BL_ASSERT(i >= 0);
        BL_ASSERT(j >= 0);
    }

    SIRec (const SIRec& rhs)
        :
        m_i(rhs.m_i),
        m_j(rhs.m_j),
        m_bx(rhs.m_bx),
        m_fbid(rhs.m_fbid) {}

    int       m_i;
    int       m_j;
    Box       m_bx;
    FillBoxId m_fbid;
};

//
// Used in caching self-intersection info for FillBoundary().
//

struct SI
{
    SI ();

    SI (const BoxArray&            ba,
        const DistributionMapping& dm,
        int                        ngrow);

    SI (const SI& rhs);

    ~SI ();

    bool operator== (const SI& rhs) const;
    bool operator!= (const SI& rhs) const;

    Array<int>          m_cache;    // Snds cached for CollectData().
    CommDataCache       m_commdata; // Yet another cache for CollectData().
    std::vector<SIRec>  m_sirec;
    BoxArray            m_ba;
    DistributionMapping m_dm;
    int                 m_ngrow;
    bool                m_reused;
};

SI::SI ()
    :
    m_ngrow(-1),
    m_reused(false)
{}

SI::SI (const BoxArray&            ba,
        const DistributionMapping& dm,
        int                        ngrow)
    :
    m_ba(ba),
    m_dm(dm),
    m_ngrow(ngrow),
    m_reused(false)
{
    BL_ASSERT(ngrow >= 0);
}

SI::SI (const SI& rhs)
    :
    m_cache(rhs.m_cache),
    m_commdata(rhs.m_commdata),
    m_sirec(rhs.m_sirec),
    m_ba(rhs.m_ba),
    m_dm(rhs.m_dm),
    m_ngrow(rhs.m_ngrow),
    m_reused(rhs.m_reused)
{}

SI::~SI () {}

bool
SI::operator== (const SI& rhs) const
{
    return m_ngrow == rhs.m_ngrow && m_ba == rhs.m_ba && m_dm == rhs.m_dm;
}

bool
SI::operator!= (const SI& rhs) const
{
    return !operator==(rhs);
}

typedef std::multimap<int,SI> SIMMap;

typedef SIMMap::iterator SIMMapIter;

static SIMMap SICache;

void
MultiFab::FlushSICache ()
{
    if (ParallelDescriptor::IOProcessor() && SICache.size())
        std::cout << "MultiFab::SICacheSize() = " << SICache.size() << std::endl;
    SICache.clear();
}

int
MultiFab::SICacheSize ()
{
    return SICache.size();
}

static
SI&
BuildFBsirec (const SI&       si,
              const MultiFab& mf)
{
    BL_ASSERT(si.m_ngrow >= 0);
    BL_ASSERT(mf.nGrow() == si.m_ngrow);
    BL_ASSERT(mf.boxArray() == si.m_ba);

    const int key = mf.nGrow() + mf.size();

    SIMMapIter it = SICache.insert(std::make_pair(key,si));

    const BoxArray&            ba     = mf.boxArray();
    const DistributionMapping& DMap   = mf.DistributionMap();
    const int                  MyProc = ParallelDescriptor::MyProc();
    std::vector<SIRec>&        sirec  = it->second.m_sirec;
    Array<int>&                cache  = it->second.m_cache;

    cache.resize(ParallelDescriptor::NProcs(),0);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const int i = mfi.index();

        std::vector< std::pair<int,Box> > isects = ba.intersections(mf[mfi].box());

        for (int ii = 0; ii < isects.size(); ii++)
        {
            const Box& bx  = isects[ii].second;
            const int  iii = isects[ii].first;

            if (i != iii)
            {
                sirec.push_back(SIRec(i,iii,bx));

                if (DMap[iii] != MyProc)
                    //
                    // If we intersect them then they'll intersect us.
                    //
                    cache[DMap[iii]] += 1;
            }
        }

        BL_ASSERT(cache[DMap[i]] == 0);
    }

    return it->second;
}

//
// Returns cached self-intersection records for MultiFab or builds them.
//

static
SI&
TheFBsirec (int             scomp,
            int             ncomp,
            const MultiFab& mf)
{
    BL_ASSERT(ncomp >  0);
    BL_ASSERT(scomp >= 0);

    static bool first             = true;
    static bool use_fb_cache      = true;
    static int  fb_cache_max_size = 25;

    if (first)
    {
        first = false;
        ParmParse pp("multifab");
        pp.query("use_fb_cache", use_fb_cache);
        pp.query("fb_cache_max_size", fb_cache_max_size);
    }

    const SI si(mf.boxArray(), mf.DistributionMap(), mf.nGrow());

    const int key = mf.nGrow() + mf.size();

    if (use_fb_cache)
    {
        std::pair<SIMMapIter,SIMMapIter> er_it = SICache.equal_range(key);
    
        for (SIMMapIter it = er_it.first; it != er_it.second; ++it)
        {
            if (it->second == si)
            {
                it->second.m_reused = true;
                //
                // Adjust the ncomp & scomp in CommData.
                //
                Array<CommData>& cd = it->second.m_commdata.theCommData();

                for (int i = 0; i < cd.size(); i++)
                {
                    cd[i].nComp(ncomp);
                    cd[i].srcComp(scomp);
                }

                return it->second;
            }
        }

        if (SICache.size() >= fb_cache_max_size)
        {
            //
            // Don't let the size of the cache get too big.
            //
            for (SIMMapIter it = SICache.begin(); it != SICache.end(); )
            {
                if (!it->second.m_reused)
                {
                    SICache.erase(it++);
                    //
                    // Only delete enough entries to stay under limit.
                    //
                    if (SICache.size() < fb_cache_max_size) break;
                }
                else
                {
                    ++it;
                }
            }

            if (SICache.size() >= fb_cache_max_size)
                //
                // Get rid of first entry which is the one with the smallest key.
                //
                SICache.erase(SICache.begin());
        }
    }
    else
    {
        SICache.clear();
    }

    return BuildFBsirec(si,mf);
}

void
MultiFab::FillBoundary (int  scomp,
                        int  ncomp,
                        bool local)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::FillBoundary(int, int)");

    if ( n_grow <= 0 ) return;

    if ( local )
    {
        //
        // Do what you can with the FABs you own.  No parallelism allowed.
        //
        const BoxArray&            ba     = boxArray();
        const DistributionMapping& DMap   = DistributionMap();
        const int                  MyProc = ParallelDescriptor::MyProc();

        for (MFIter mfi(*this); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();

            std::vector< std::pair<int,Box> > isects = ba.intersections((*this)[mfi].box());

            for (int ii = 0; ii < isects.size(); ii++)
            {
                const Box& bx  = isects[ii].second;
                const int  iii = isects[ii].first;

                if (i != iii && DMap[iii] == MyProc)
                {
                    (*this)[i].copy((*this)[iii], bx, scomp, bx, scomp, ncomp);
                }
            }
        }
    }
    else
    {
        MultiFabCopyDescriptor mfcd;
        SI&                    si   = TheFBsirec(scomp,ncomp,*this);
        const MultiFabId       mfid = mfcd.RegisterMultiFab(this);
        //
        // Add boxes we need to collect.
        //
        for (unsigned int i = 0; i < si.m_sirec.size(); i++)
        {
            si.m_sirec[i].m_fbid = mfcd.AddBox(mfid,
                                               si.m_sirec[i].m_bx,
                                               0,
                                               si.m_sirec[i].m_j,
                                               scomp,
                                               scomp,
                                               ncomp);
        }

        mfcd.CollectData(&si.m_cache,&si.m_commdata);

        for (unsigned int i = 0; i < si.m_sirec.size(); i++)
        {
            BL_ASSERT(DistributionMap()[si.m_sirec[i].m_i] == ParallelDescriptor::MyProc());
            //
            // Directly fill the FAB.
            //
            mfcd.FillFab(mfid,si.m_sirec[i].m_fbid,(*this)[si.m_sirec[i].m_i]);
        }
    }

}
