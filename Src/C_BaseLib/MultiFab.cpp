//
// $Id: MultiFab.cpp,v 1.52 2001-07-17 23:02:25 lijewski Exp $
//

#include <algorithm>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <list>

#include <BLassert.H>
#include <MultiFab.H>
#include <ParallelDescriptor.H>

#ifdef BL3_PROFILING
#include <BoxLib3/Profiler.H>
#endif

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

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

    for (MultiFabIterator mfi(dst); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi,src);

        BL_ASSERT(mfi.validbox() == dmfi.validbox());

#ifndef BL_NAMESPACE
        Box bx = ::grow(mfi.validbox(),nghost);
#else
        Box bx = BL_NAMESPACE::grow(mfi.validbox(),nghost);
#endif

        if (bx.ok())
            mfi().copy(dmfi(), bx, srccomp, bx, dstcomp, numcomp);
    }
}

void
MultiFab::FillBoundary ()
{
    FillBoundary(0, n_comp);
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

MultiFab::MultiFab () {}

MultiFab::MultiFab (const BoxArray& bxs,
                    int             ncomp,
                    int             ngrow,
                    FabAlloc        alloc)
    :
    FabArray<Real,FArrayBox>(bxs,ncomp,ngrow,alloc)
{}

//
// This isn't inlined as it's virtual.
//
MultiFab::~MultiFab () {}

void
MultiFab::probe (std::ostream& os,
                 IntVect&      pt)
{
    Real  dat[20];
    int prec = os.precision(14);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        if (mfi.validbox().contains(pt))
        {
            mfi().getVal(dat,pt);

            os << "point "
               << pt
               << " in box "
               << mfi.validbox()
               << " data = ";
            for (int i = 0, N = mfi().nComp(); i < N; i++)
                os << ' ' << std::setw(20) << dat[i];
            os << '\n';
        }
    }
    os.precision(prec);

    if (os.fail())
        BoxLib::Error("MultiFab::probe(ostream&,IntVect&) failed");
}

Real
MultiFab::min (int comp,
               int nghost) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);

#ifdef  BL_USE_FLOAT
    Real mn = FLT_MAX;
#else
    Real mn = DBL_MAX;
#endif

    for (ConstMultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
#ifndef BL_NAMESPACE
        mn = std::min(mn,mfi().min(::grow(mfi.validbox(),nghost),comp));
#else
        mn = std::min(mn,mfi().min(BL_NAMESPACE::grow(mfi.validbox(),nghost),comp));
#endif
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

#ifdef  BL_USE_FLOAT
    Real mn = FLT_MAX;
#else
    Real mn = DBL_MAX;
#endif

    for (ConstMultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
#ifndef BL_NAMESPACE
        Box b = ::grow(mfi.validbox(),nghost) & region;
#else
        Box b = BL_NAMESPACE::grow(mfi.validbox(),nghost) & region;
#endif

        if (b.ok())
            mn = std::min(mn,mfi().min(b,comp));
    }

    ParallelDescriptor::ReduceRealMin(mn);

    return mn;
}

Real
MultiFab::max (int comp,
               int nghost) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);

#ifdef  BL_USE_FLOAT
    Real mn = -FLT_MAX;
#else
    Real mn = -DBL_MAX;
#endif

    for (ConstMultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
#ifndef BL_NAMESPACE
        mn = std::max(mn,mfi().max(::grow(mfi.validbox(),nghost),comp));
#else
        mn = Max(mn,mfi().max(BL_NAMESPACE::grow(mfi.validbox(),nghost),comp));
#endif
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

#ifdef  BL_USE_FLOAT
    Real mn = -FLT_MAX;
#else
    Real mn = -DBL_MAX;
#endif

    for (ConstMultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
#ifndef BL_NAMESPACE
        Box b = ::grow(mfi.validbox(),nghost) & region;
#else
        Box b = BL_NAMESPACE::grow(mfi.validbox(),nghost) & region;
#endif

        if (b.ok())
            mn = std::max(mn,mfi().max(b,comp));
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
#ifndef NDEBUG
    int lst_comp = strt_comp + num_comp - 1;
#endif
    BL_ASSERT(lst_comp < n_comp && lst_comp < mf.n_comp);
    BL_ASSERT(nghost <= n_grow && nghost <= mf.n_grow);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi, mf);

#ifndef BL_NAMESPACE
        Box bx = ::grow(mfi.validbox(),nghost);
#else
        Box bx = BL_NAMESPACE::grow(mfi.validbox(),nghost);
#endif

        mfi().minus(dmfi(), bx, strt_comp, strt_comp, num_comp);
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

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
#ifndef BL_NAMESPACE
        mfi().plus(val,::grow(mfi.validbox(),nghost),comp,num_comp);
#else
        mfi().plus(val,BL_NAMESPACE::grow(mfi.validbox(),nghost),comp,num_comp);
#endif
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

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
#ifndef BL_NAMESPACE
        Box b = ::grow(mfi.validbox(),nghost) & region;
#else
        Box b = BL_NAMESPACE::grow(mfi.validbox(),nghost) & region;
#endif

        if (b.ok())
            mfi().plus(val,b,comp,num_comp);
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
#ifndef NDEBUG
    int lst_comp = strt_comp + num_comp - 1;
#endif
    BL_ASSERT(lst_comp < n_comp && lst_comp < mf.n_comp);
    BL_ASSERT(nghost <= n_grow && nghost <= mf.n_grow);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi, mf);

#ifndef BL_NAMESPACE
        Box bx = ::grow(mfi.validbox(),nghost);
#else
        Box bx = BL_NAMESPACE::grow(mfi.validbox(),nghost);
#endif

        mfi().plus(dmfi(), bx, strt_comp, strt_comp, num_comp);
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

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
#ifndef BL_NAMESPACE
        mfi().mult(val, ::grow(mfi.validbox(),nghost),comp,num_comp);
#else
        mfi().mult(val,BL_NAMESPACE::grow(mfi.validbox(),nghost),comp,num_comp);
#endif
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

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
#ifndef BL_NAMESPACE
        Box b = ::grow(mfi.validbox(),nghost) & region;
#else
        Box b = BL_NAMESPACE::grow(mfi.validbox(),nghost) & region;
#endif

        if (b.ok())
            mfi().mult(val, b, comp, num_comp);
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

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
#ifndef BL_NAMESPACE
        mfi().invert(numerator, ::grow(mfi.validbox(),nghost),comp,num_comp);
#else
        mfi().invert(numerator,BL_NAMESPACE::grow(mfi.validbox(),nghost),comp,num_comp);
#endif
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

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
#ifndef BL_NAMESPACE
        Box b = ::grow(mfi.validbox(),nghost) & region;
#else
        Box b = BL_NAMESPACE::grow(mfi.validbox(),nghost) & region;
#endif

        if (b.ok())
            mfi().invert(numerator,b,comp,num_comp);
    }
}

void
MultiFab::negate (int comp,
                  int num_comp,
                  int nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
#ifndef BL_NAMESPACE
        mfi().negate(::grow(mfi.validbox(),nghost),comp,num_comp);
#else
        mfi().negate(BL_NAMESPACE::grow(mfi.validbox(),nghost),comp,num_comp);
#endif
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

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
#ifndef BL_NAMESPACE
        Box b = ::grow(mfi.validbox(),nghost) & region;
#else
        Box b = BL_NAMESPACE::grow(mfi.validbox(),nghost) & region;
#endif

        if (b.ok())
            mfi().negate(b,comp,num_comp);
    }
}

void
linInterpAddBox (MultiFabCopyDescriptor& fabCopyDesc,
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
linInterpFillFab (MultiFabCopyDescriptor& fabCopyDesc,
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
        dest1.setVal(1.e30);
        FArrayBox dest2(dest.box(), dest.nComp());
        dest2.setVal(1.e30);
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

    SI (const BoxArray& ba,
        int             scomp,
        int             ncomp,
        int             ngrow);

    SI (const SI& rhs);

    ~SI ();

    bool operator== (const SI& rhs) const;
    bool operator!= (const SI& rhs) const;

    Array<int>    m_cache;    // Snds cached for CollectData().
    CommDataCache m_commdata; // Yet another cache for CollectData().
    std::vector<SIRec> m_sirec;
    BoxArray      m_ba;
    int           m_scomp;
    int           m_ncomp;
    int           m_ngrow;
};

inline
SI::SI ()
    :
    m_scomp(-1),
    m_ncomp(-1),
    m_ngrow(-1)
{}

inline
SI::SI (const BoxArray& ba,
        int             scomp,
        int             ncomp,
        int             ngrow)
    :
    m_ba(ba),
    m_scomp(scomp),
    m_ncomp(ncomp),
    m_ngrow(ngrow)
{
    BL_ASSERT(ncomp >  0);
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ngrow >= 0);
}

inline
SI::SI (const SI& rhs)
    :
    m_cache(rhs.m_cache),
    m_commdata(rhs.m_commdata),
    m_sirec(rhs.m_sirec),
    m_ba(rhs.m_ba),
    m_scomp(rhs.m_scomp),
    m_ncomp(rhs.m_ncomp),
    m_ngrow(rhs.m_ngrow)
{}

inline
SI::~SI () {}

inline
bool
SI::operator== (const SI& rhs) const
{
    return
        m_scomp == rhs.m_scomp &&
        m_ncomp == rhs.m_ncomp &&
        m_ngrow == rhs.m_ngrow &&
        m_ba    == rhs.m_ba;
}

inline
bool
SI::operator!= (const SI& rhs) const
{
    return !operator==(rhs);
}

//
// A useful typedef.
//
typedef std::list<SI> SIList;

//
// Cache of SI info.
//
static SIList SICache;

void
MultiFab::FlushSICache ()
{
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
    BL_ASSERT(si.m_ncomp >  0);
    BL_ASSERT(si.m_scomp >= 0);
    BL_ASSERT(si.m_ngrow >= 0);
    BL_ASSERT(mf.nGrow() == si.m_ngrow);
    BL_ASSERT(mf.boxArray() == si.m_ba);
    //
    // Insert new ones at beginning of list.
    //
    SICache.push_front(si);

    //cout << "*** FB Cache Size = " << SICache.size() << endl;

    const BoxArray&            ba     = mf.boxArray();
    const DistributionMapping& DMap   = mf.DistributionMap();
    const int                  MyProc = ParallelDescriptor::MyProc();
    std::vector<SIRec>&        sirec  = SICache.front().m_sirec;
    Array<int>&                cache  = SICache.front().m_cache;

    cache.resize(ParallelDescriptor::NProcs(),0);

    for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
    {
        const int i = mfi.index();

        for (int j = 0; j < mf.length(); j++)
        {
            if (i != j)
            {
                if (ba[j].intersects(mfi().box()))
                {
                    Box bx = ba[j] & mfi().box();

                    sirec.push_back(SIRec(i,j,bx));

                    if (DMap[j] != MyProc)
                        //
                        // If we intersect them then they'll intersect us.
                        //
                        cache[DMap[j]] += 1;
                }
            }
        }

        BL_ASSERT(cache[DMap[i]] == 0);
    }

    return SICache.front();
}

//
// Returns cached self-intersection records for MultiFab or builds them.
//

inline
SI&
TheFBsirec (int             scomp,
            int             ncomp,
            const MultiFab& mf)
{
    BL_ASSERT(ncomp >  0);
    BL_ASSERT(scomp >= 0);

    const SI si(mf.boxArray(), scomp, ncomp, mf.nGrow());
    
    for (SIList::iterator it = SICache.begin(); it != SICache.end(); ++it)
        if (*it == si)
            return *it;

    return BuildFBsirec(si,mf);
}

void
MultiFab::FillBoundary (int scomp,
                        int ncomp)
{
#ifdef BL3_PROFILING
  BL3_PROFILE(BL3_PROFILE_THIS_NAME() + "::FillBoundary(int, int)");
#endif
    static RunStats stats("fill_boundary");

    stats.start();

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

    stats.end();
}

#ifdef BL_NAMESPACE
}
#endif
