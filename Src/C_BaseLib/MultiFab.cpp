//
// $Id: MultiFab.cpp,v 1.90 2008-10-29 16:50:44 sepp Exp $
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
MultiFab::Multiply (MultiFab&       dst,
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
            dst[mfi].mult(src[mfi], bx, bx, srccomp, dstcomp, numcomp);
    }
}

void
MultiFab::Divide (MultiFab&       dst,
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
            dst[mfi].divide(src[mfi], bx, bx, srccomp, dstcomp, numcomp);
    }
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

Real
MultiFab::norm0 (int comp) const
{
    Real nm0 = -std::numeric_limits<Real>::max();

    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        nm0 = std::max(nm0, get(mfi).norm(mfi.validbox(), 0, comp, 1));
    }

    ParallelDescriptor::ReduceRealMax(nm0);

    return nm0;
}

Real
MultiFab::norm2 (int comp) const
{
    Real nm2 = 0.0;
    Real nm_grid;

    for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
        nm_grid = get(mfi).norm(mfi.validbox(), 2, comp, 1);
        nm2 += nm_grid*nm_grid;
    }

    ParallelDescriptor::ReduceRealSum(nm2);

    nm2 = std::sqrt(nm2);

    return nm2;
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

void
MultiFab::FillBoundary (int  scomp,
                        int  ncomp,
                        bool local)
{
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
        FabArray<FArrayBox>::FillBoundary(scomp,ncomp);
    }
}

void
MultiFab::FillBoundary (bool local)
{
    FillBoundary(0, n_comp, local);
}

