//BL_COPYRIGHT_NOTICE

//
// $Id: MultiFab.cpp,v 1.40 1999-05-10 17:18:46 car Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <list>
using std::list;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;
#else
#include <float.h>
#include <iostream.h>
#include <iomanip.h>
#include <list.h>
#endif

#include <BLassert.H>
#include <Misc.H>
#include <MultiFab.H>
#include <ParallelDescriptor.H>

void
MultiFab::Copy (MultiFab&       dst,
                const MultiFab& src,
                int             srccomp,
                int             dstcomp,
                int             numcomp,
                int             nghost)
{
    BLassert(dst.boxArray() == src.boxArray());
    BLassert(dst.distributionMap == src.distributionMap);
    BLassert(dst.nGrow() >= nghost && src.nGrow() >= nghost);

    for (MultiFabIterator mfi(dst); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi,src);

        BLassert(mfi.validbox() == dmfi.validbox());

        Box bx = ::grow(mfi.validbox(),nghost);

        if (bx.ok())
            mfi().copy(dmfi(), bx, srccomp, bx, dstcomp, numcomp);
    }
}

MultiFab::MultiFab ()
    :
    m_FB_mfcd(0),
    m_FB_scomp(-1),
    m_FB_ncomp(-1),
    m_FPB_mfcd(0),
    m_FPB_scomp(-1),
    m_FPB_ncomp(-1),
    m_FPB_noovlp(false)
{}

MultiFab::MultiFab (const BoxArray& bxs,
                    int             ncomp,
                    int             ngrow,
                    FabAlloc        alloc)
    :
    FabArray<Real,FArrayBox>(bxs,ncomp,ngrow,alloc),
    m_FB_mfcd(0),
    m_FB_scomp(-1),
    m_FB_ncomp(-1),
    m_FPB_mfcd(0),
    m_FPB_scomp(-1),
    m_FPB_ncomp(-1),
    m_FPB_noovlp(false)
{}

//
// This isn't inlined as it's virtual.
//
MultiFab::~MultiFab ()
{
    delete m_FB_mfcd;
    delete m_FPB_mfcd;
}

void
MultiFab::probe (ostream& os,
                 IntVect& pt)
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
                os << ' ' << setw(20) << dat[i];
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
    BLassert(nghost >= 0 && nghost <= n_grow);

#ifdef BL_USE_DOUBLE
    Real mn = FLT_MAX;
#elif  BL_USE_FLOAT
    Real mn = DBL_MAX;
#endif

    for (ConstMultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        mn = Min(mn,mfi().min(::grow(mfi.validbox(),nghost),comp));
    }

    ParallelDescriptor::ReduceRealMin(mn);

    return mn;
}

Real
MultiFab::min (const Box& region,
               int        comp,
               int        nghost) const
{
    BLassert(nghost >= 0 && nghost <= n_grow);

#ifdef BL_USE_DOUBLE
    Real mn = FLT_MAX;
#elif  BL_USE_FLOAT
    Real mn = DBL_MAX;
#endif

    for (ConstMultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b = ::grow(mfi.validbox(),nghost) & region;

        if (b.ok())
            mn = Min(mn,mfi().min(b,comp));
    }

    ParallelDescriptor::ReduceRealMin(mn);

    return mn;
}

Real
MultiFab::max (int comp,
               int nghost) const
{
    BLassert(nghost >= 0 && nghost <= n_grow);

#ifdef BL_USE_DOUBLE
    Real mn = FLT_MIN;
#elif  BL_USE_FLOAT
    Real mn = DBL_MIN;
#endif

    for (ConstMultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        mn = Max(mn,mfi().max(::grow(mfi.validbox(),nghost),comp));
    }

    ParallelDescriptor::ReduceRealMax(mn);

    return mn;
}

Real
MultiFab::max (const Box& region,
               int        comp,
               int        nghost) const
{
    BLassert(nghost >= 0 && nghost <= n_grow);

#ifdef BL_USE_DOUBLE
    Real mn = FLT_MIN;
#elif  BL_USE_FLOAT
    Real mn = DBL_MIN;
#endif

    for (ConstMultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b = ::grow(mfi.validbox(),nghost) & region;

        if (b.ok())
            mn = Max(mn,mfi().max(b,comp));
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
    BLassert(boxarray == mf.boxarray);
    BLassert(strt_comp >= 0);
#ifndef NDEBUG
    int lst_comp = strt_comp + num_comp - 1;
#endif
    BLassert(lst_comp < n_comp && lst_comp < mf.n_comp);
    BLassert(nghost <= n_grow && nghost <= mf.n_grow);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi, mf);

        Box bx = ::grow(mfi.validbox(),nghost);

        mfi().minus(dmfi(), bx, strt_comp, strt_comp, num_comp);
    }
}

void
MultiFab::plus (Real val,
                int  comp,
                int  num_comp,
                int  nghost)
{
    BLassert(nghost >= 0 && nghost <= n_grow);
    BLassert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        mfi().plus(val,::grow(mfi.validbox(),nghost),comp,num_comp);
    }
}

void
MultiFab::plus (Real       val,
                const Box& region,
                int        comp,
                int        num_comp,
                int        nghost)
{
    BLassert(nghost >= 0 && nghost <= n_grow);
    BLassert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b = ::grow(mfi.validbox(),nghost) & region;

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
    BLassert(boxarray == mf.boxarray);
    BLassert(strt_comp >= 0);
#ifndef NDEBUG
    int lst_comp = strt_comp + num_comp - 1;
#endif
    BLassert(lst_comp < n_comp && lst_comp < mf.n_comp);
    BLassert(nghost <= n_grow && nghost <= mf.n_grow);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi, mf);

        Box bx = ::grow(mfi.validbox(),nghost);

        mfi().plus(dmfi(), bx, strt_comp, strt_comp, num_comp);
    }
}

void
MultiFab::mult (Real val,
                int  comp,
                int  num_comp,
                int  nghost)
{
    BLassert(nghost >= 0 && nghost <= n_grow);
    BLassert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        mfi().mult(val,::grow(mfi.validbox(),nghost),comp,num_comp);
    }
}

void
MultiFab::mult (Real       val,
                const Box& region,
                int        comp,
                int        num_comp,
                int        nghost)
{
    BLassert(nghost >= 0 && nghost <= n_grow);
    BLassert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b = ::grow(mfi.validbox(),nghost) & region;

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
    BLassert(nghost >= 0 && nghost <= n_grow);
    BLassert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        mfi().invert(numerator,::grow(mfi.validbox(),nghost),comp,num_comp);
    }
}

void
MultiFab::invert (Real       numerator,
                  const Box& region,
                  int        comp,
                  int        num_comp,
                  int        nghost)
{
    BLassert(nghost >= 0 && nghost <= n_grow);
    BLassert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b = ::grow(mfi.validbox(),nghost) & region;

        if (b.ok())
            mfi().invert(numerator,b,comp,num_comp);
    }
}

void
MultiFab::negate (int comp,
                  int num_comp,
                  int nghost)
{
    BLassert(nghost >= 0 && nghost <= n_grow);
    BLassert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        mfi().negate(::grow(mfi.validbox(),nghost),comp,num_comp);
    }
}

void
MultiFab::negate (const Box& region,
                  int        comp,
                  int        num_comp,
                  int        nghost)
{
    BLassert(nghost >= 0 && nghost <= n_grow);
    BLassert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b = ::grow(mfi.validbox(),nghost) & region;

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

    BLassert(t>t1-teps && (extrap || t < t2+teps));

    if (t < t1+teps)
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

    BLassert(t>t1-teps && (extrap || t < t2+teps));

    if (t < t1+teps)
    {
        fabCopyDesc.FillFab(faid1, fillBoxIds[0], dest);
    }
    else if (t > t2-teps && t < t2+teps)
    {
        fabCopyDesc.FillFab(faid2, fillBoxIds[0], dest);
    }
    else
    {
        BLassert(dest_comp + num_comp <= dest.nComp());

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
        BLassert(i >= 0);
        BLassert(j >= 0);
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
    bool operator!= (const SI& rhs) const { return !operator==(rhs); }

    Array<int>    m_cache;    // Snds cached for CollectData().
    CommDataCache m_commdata; // Yet another cache for CollectData().
    vector<SIRec> m_sirec;
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
    BLassert(ncomp >  0);
    BLassert(scomp >= 0);
    BLassert(ngrow >= 0);
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

//
// A useful typedef.
//
typedef list<SI> SIList;

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
    BLassert(si.m_ncomp >  0);
    BLassert(si.m_scomp >= 0);
    BLassert(si.m_ngrow >= 0);
    BLassert(mf.nGrow() == si.m_ngrow);
    BLassert(mf.boxArray() == si.m_ba);
    //
    // Insert new ones at beginning of list.
    //
    SICache.push_front(si);

    //cout << "*** FB Cache Size = " << SICache.size() << endl;

    const BoxArray&            ba     = mf.boxArray();
    const DistributionMapping& DMap   = mf.DistributionMap();
    const int                  MyProc = ParallelDescriptor::MyProc();
    vector<SIRec>&             sirec  = SICache.front().m_sirec;
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

        BLassert(cache[DMap[i]] == 0);
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
    BLassert(ncomp >  0);
    BLassert(scomp >= 0);

    const SI si(mf.boxArray(), scomp, ncomp, mf.nGrow());
    
    for (SIList::iterator it = SICache.begin(); it != SICache.end(); ++it)
    {
        if (*it == si)
        {
            return *it;
        }
    }

    return BuildFBsirec(si,mf);
}

void
MultiFab::FillBoundary (int scomp,
                        int ncomp)
{
    static RunStats stats("fill_boundary");

    stats.start();

    MultiFabCopyDescriptor& mfcd = theFBmfcd(scomp,ncomp);
    const MultiFabId        mfid = 0;
    SI&                     si   = TheFBsirec(scomp,ncomp,*this);

    if (mfcd.nFabComTags() == 0)
    {
        //
        // Add boxes we need to collect, if we haven't already done so.
        //
        for (int i = 0; i < si.m_sirec.size(); i++)
        {
            si.m_sirec[i].m_fbid = mfcd.AddBox(mfid,
                                               si.m_sirec[i].m_bx,
                                               0,
                                               si.m_sirec[i].m_j,
                                               scomp,
                                               scomp,
                                               ncomp);
        }
    }

    mfcd.CollectData(&si.m_cache,&si.m_commdata);

    for (int i = 0; i < si.m_sirec.size(); i++)
    {
        BLassert(DistributionMap()[si.m_sirec[i].m_i] == ParallelDescriptor::MyProc());
        //
        // Directly fill the FAB.
        //
        mfcd.FillFab(mfid,si.m_sirec[i].m_fbid,(*this)[si.m_sirec[i].m_i]);
    }

    stats.end();
}
