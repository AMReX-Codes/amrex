
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
#include <BLProfiler.H>
#include <ParmParse.H>
#include <PArray.H>

//
// Set default values in Initialize()!!!
//
bool MultiFab::check_for_nan;
bool MultiFab::check_for_inf;

namespace
{
    bool initialized = false;
}

MultiFabCopyDescriptor::MultiFabCopyDescriptor ()
    :
    FabArrayCopyDescriptor<FArrayBox>() {}

MultiFabCopyDescriptor::~MultiFabCopyDescriptor () {}

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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(dst,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);

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
// don't have to    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrow() >= nghost); // && src.nGrow() >= nghost);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(dst,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);

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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(dst,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);

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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(dst,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);

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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(dst,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);

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
MultiFab::Initialize ()
{
    if (initialized) return;
    //
    // Set initial values!!!
    //
    MultiFab::check_for_nan = false;
    MultiFab::check_for_inf = false;

    ParmParse pp("multifab");

    pp.query("check_for_nan", check_for_nan);
    pp.query("check_for_inf", check_for_inf);

    BoxLib::ExecOnFinalize(MultiFab::Finalize);

    initialized = true;
}

void
MultiFab::Finalize ()
{
    initialized = false;
}

MultiFab::MultiFab ()
{
    Initialize();
}

MultiFab::MultiFab (const BoxArray& bxs,
                    int             ncomp,
                    int             ngrow,
                    FabAlloc        alloc,
                    const IntVect&  nodal)
    :
    FabArray<FArrayBox>(bxs,ncomp,ngrow,alloc,nodal)
{
    Initialize();

    if ((check_for_nan || check_for_inf) && alloc == Fab_allocate) setVal(0);
}

MultiFab::MultiFab (const BoxArray&            bxs,
                    int                        ncomp,
                    int                        ngrow,
                    const DistributionMapping& dm,
                    FabAlloc                   alloc,
                    const IntVect&             nodal)
    :
    FabArray<FArrayBox>(bxs,ncomp,ngrow,dm,alloc,nodal)
{
    Initialize();

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
                  FabAlloc        alloc,
		  const IntVect&  nodal)
{
    this->FabArray<FArrayBox>::define(bxs,nvar,ngrow,alloc,nodal);

    if ((check_for_nan || check_for_inf) && alloc == Fab_allocate) setVal(0);
}

void
MultiFab::define (const BoxArray&            bxs,
                  int                        nvar,
                  int                        ngrow,
                  const DistributionMapping& dm,
                  FabAlloc                   alloc,
		  const IntVect&             nodal)
{
    this->FabArray<FArrayBox>::define(bxs,nvar,ngrow,dm,alloc,nodal);

    if ((check_for_nan || check_for_inf) && alloc == Fab_allocate) setVal(0);
}

bool 
MultiFab::contains_nan (int scomp,
                        int ncomp,
                        int ngrow,
			bool local) const
{
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(scomp + ncomp <= nComp());
    BL_ASSERT(ncomp >  0 && ncomp <= nComp());
    BL_ASSERT(ngrow >= 0 && ngrow <= nGrow());

    bool r = false;

#ifdef _OPENMP
#pragma omp parallel reduction(|:r)
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.growntilebox(ngrow);
	
	if (this->FabArray<FArrayBox>::get(mfi).contains_nan(bx,scomp,ncomp))
	    r = true;
    }

    if (!local)
	ParallelDescriptor::ReduceBoolOr(r);

    return r;
}

bool 
MultiFab::contains_nan (bool local) const
{
    return contains_nan(0,nComp(),nGrow(),local);
}

bool 
MultiFab::contains_inf (int scomp,
                        int ncomp,
                        int ngrow,
			bool local) const
{
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(scomp + ncomp <= nComp());
    BL_ASSERT(ncomp >  0 && ncomp <= nComp());
    BL_ASSERT(ngrow >= 0 && ngrow <= nGrow());

    bool r = false;

#ifdef _OPENMP
#pragma omp parallel reduction(|:r)
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.growntilebox(ngrow);
	
	if (this->FabArray<FArrayBox>::get(mfi).contains_inf(bx,scomp,ncomp))
	    r = true;
    }

    if (!local)
	ParallelDescriptor::ReduceBoolOr(r);

    return r;
}

bool 
MultiFab::contains_inf (bool local) const
{
    return contains_inf(0,nComp(),nGrow(),local);
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
MultiFab::operator[] (int K) const
{
    BL_ASSERT(defined(K));

    const FArrayBox& fab = this->FabArray<FArrayBox>::get(K);

    if (check_for_nan && fab.contains_nan())
        AbortOnNaN(fab);

    if (check_for_inf && fab.contains_inf())
        AbortOnInf(fab);

    return fab;
}

FArrayBox&
MultiFab::operator[] (int K)
{
    BL_ASSERT(defined(K));

    FArrayBox& fab = this->FabArray<FArrayBox>::get(K);

    if (check_for_nan && fab.contains_nan())
        AbortOnNaN(fab);

    if (check_for_inf && fab.contains_inf())
        AbortOnInf(fab);

    return fab;
}

Real
MultiFab::min (int comp,
               int nghost,
	       bool local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);

    Real mn = std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel reduction(min:mn)
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.growntilebox(nghost);
	mn = std::min(mn, get(mfi).min(bx,comp));
    }

    if (!local)
	ParallelDescriptor::ReduceRealMin(mn);

    return mn;
}

Real
MultiFab::min (const Box& region,
               int        comp,
               int        nghost,
	       bool       local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);

    Real mn = std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel reduction(min:mn)
#endif
    for ( MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
	const Box& b = mfi.growntilebox(nghost) & region;
	
	if (b.ok())
	    mn = std::min(mn, get(mfi).min(b,comp));
    }

    if (!local)
	ParallelDescriptor::ReduceRealMin(mn);

    return mn;
}

Real
MultiFab::max (int comp,
               int nghost,
	       bool local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);

    Real mx = -std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel reduction(max:mx)
#endif
    for ( MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.growntilebox(nghost);
	mx = std::max(mx, get(mfi).max(bx,comp));
    }

    if (!local)
	ParallelDescriptor::ReduceRealMax(mx);

    return mx;
}

Real
MultiFab::max (const Box& region,
               int        comp,
               int        nghost,
	       bool       local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);

    Real mx = -std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel reduction(max:mx)
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
	const Box& b = mfi.growntilebox(nghost) & region;

	if (b.ok())
	    mx = std::max(mx, get(mfi).max(b,comp));
    }
	
    if (!local)
	ParallelDescriptor::ReduceRealMax(mx);

    return mx;
}

IntVect
MultiFab::minIndex (int comp,
                    int nghost) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);

    IntVect loc;

    Real mn = std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	Real priv_mn = std::numeric_limits<Real>::max();
	IntVect priv_loc;

	for (MFIter mfi(*this); mfi.isValid(); ++mfi)
	{
	    const Box& box = BoxLib::grow(mfi.validbox(),nghost);
	    const Real lmn = get(mfi).min(box,comp);

	    if (lmn < priv_mn)
	    {
		priv_mn  = lmn;
		priv_loc = get(mfi).minIndex(box,comp);
	    }
	}
#ifdef _OPENMP
#pragma omp critical (multifab_minindex)
#endif
	{
	    if (priv_mn < mn) {
		mn = priv_mn;
		loc = priv_loc;
	    }
	}
    }

    const int NProcs = ParallelDescriptor::NProcs();

    if (NProcs > 1)
    {
        Array<Real> mns(1);
        Array<int>  locs(1);

        if (ParallelDescriptor::IOProcessor())
        {
            mns.resize(NProcs);
            locs.resize(NProcs*BL_SPACEDIM);
        }

        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::Gather(&mn, 1, mns.dataPtr(), 1, IOProc);

        BL_ASSERT(sizeof(IntVect) == sizeof(int)*BL_SPACEDIM);

        ParallelDescriptor::Gather(loc.getVect(), BL_SPACEDIM, locs.dataPtr(), BL_SPACEDIM, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            mn  = mns[0];
            loc = IntVect(D_DECL(locs[0],locs[1],locs[2]));

            for (int i = 1; i < NProcs; i++)
            {
                if (mns[i] < mn)
                {
                    mn = mns[i];

                    const int j = BL_SPACEDIM * i;

                    loc = IntVect(D_DECL(locs[j+0],locs[j+1],locs[j+2]));
                }
            }
        }

        ParallelDescriptor::Bcast(const_cast<int*>(loc.getVect()), BL_SPACEDIM, IOProc);
    }

    return loc;
}

IntVect
MultiFab::maxIndex (int comp,
                    int nghost) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);

    IntVect loc;

    Real mx = -std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	Real priv_mx = -std::numeric_limits<Real>::max();
	IntVect priv_loc;
	
	for (MFIter mfi(*this); mfi.isValid(); ++mfi)
	{
	    const Box& box = BoxLib::grow(mfi.validbox(),nghost);
	    const Real lmx = get(mfi).max(box,comp);

	    if (lmx > priv_mx)
	    {
		priv_mx  = lmx;
		priv_loc = get(mfi).maxIndex(box,comp);
	    }
	}
#ifdef _OPENMP
#pragma omp critical (multifab_maxindex)
#endif
	{
	    if (priv_mx > mx) {
		mx = priv_mx;
		loc = priv_loc;
	    }
	}
    }

    const int NProcs = ParallelDescriptor::NProcs();

    if (NProcs > 1)
    {
        Array<Real> mxs(1);
        Array<int>  locs(1);

        if (ParallelDescriptor::IOProcessor())
        {
            mxs.resize(NProcs);
            locs.resize(NProcs*BL_SPACEDIM);
        }

        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::Gather(&mx, 1, mxs.dataPtr(), 1, IOProc);

        BL_ASSERT(sizeof(IntVect) == sizeof(int)*BL_SPACEDIM);

        ParallelDescriptor::Gather(loc.getVect(), BL_SPACEDIM, locs.dataPtr(), BL_SPACEDIM, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            mx  = mxs[0];
            loc = IntVect(D_DECL(locs[0],locs[1],locs[2]));

            for (int i = 1; i < NProcs; i++)
            {
                if (mxs[i] > mx)
                {
                    mx = mxs[i];

                    const int j = BL_SPACEDIM * i;

                    loc = IntVect(D_DECL(locs[j+0],locs[j+1],locs[j+2]));
                }
            }
        }

        ParallelDescriptor::Bcast(const_cast<int*>(loc.getVect()), BL_SPACEDIM, IOProc);
    }

    return loc;
}

Real
MultiFab::norm0 (int comp, const BoxArray& ba, bool local) const
{
    Real nm0 = -std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel reduction(max:nm0)
#endif
    {
	std::vector< std::pair<int,Box> > isects;

	for (MFIter mfi(*this); mfi.isValid(); ++mfi)
	{
	    ba.intersections(mfi.validbox(),isects);

	    for (int i = 0, N = isects.size(); i < N; i++)
	    {
		nm0 = std::max(nm0, get(mfi).norm(isects[i].second, 0, comp, 1));
	    }
	}
    }
 
    if (!local)
	ParallelDescriptor::ReduceRealMax(nm0);
 
    return nm0;
}

Real
MultiFab::norm0 (int comp, bool local) const
{
    Real nm0 = -std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel reduction(max:nm0)
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
	nm0 = std::max(nm0, get(mfi).norm(mfi.tilebox(), 0, comp, 1));
    }

    if (!local)
	ParallelDescriptor::ReduceRealMax(nm0);

    return nm0;
}

Array<Real>
MultiFab::norm0 (const Array<int>& comps, bool local) const
{
    int n = comps.size();
    const Real rmax = std::numeric_limits<Real>::max();
    Array<Real> nm0(n, -rmax);

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif
    PArray< Array<Real> > priv_nm0(nthreads, PArrayManage);
    for (int i=0; i<nthreads; i++) {
	priv_nm0.set(i, new Array<Real>(n, -rmax));
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
	int tid = omp_get_thread_num();
#else
	int tid = 0;
#endif
	for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
	{
            for (int i=0; i<n; i++) {
	        priv_nm0[tid][i] = std::max(priv_nm0[tid][i], get(mfi).norm(mfi.tilebox(), 0, comps[i], 1));
            }
        }
#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
#endif
	for (int i=0; i<n; i++) {
            for (int it=0; it<nthreads; it++) {
	        nm0[i] = std::max(priv_nm0[it][i], nm0[i]);
            }	    
	}
    }

    if (!local)
	ParallelDescriptor::ReduceRealMax(nm0.dataPtr(), n);

    return nm0;
}

Real
MultiFab::norm2 (int comp) const
{
    Real nm2 = 0.e0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:nm2)
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        const Real nm_grid = get(mfi).norm(mfi.tilebox(), 2, comp, 1);

        nm2 += nm_grid*nm_grid;
    }

    ParallelDescriptor::ReduceRealSum(nm2);

    nm2 = std::sqrt(nm2);

    return nm2;
}

Array<Real>
MultiFab::norm2 (const Array<int>& comps) const
{
    int n = comps.size();
    Array<Real> nm2(n, 0.e0);

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif
    PArray< Array<Real> > priv_nm2(nthreads, PArrayManage);
    for (int i=0; i<nthreads; i++) {
	priv_nm2.set(i, new Array<Real>(n, 0.0));
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
	int tid = omp_get_thread_num();
#else
	int tid = 0;
#endif
	for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
	{
            for (int i=0; i<n; i++) {
	        const Real nm_grid = get(mfi).norm(mfi.tilebox(), 2, comps[i], 1);
		priv_nm2[tid][i] += nm_grid*nm_grid;
            }
        }
#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
#endif
	for (int i=0; i<n; i++) {
	    for (int it=0; it<nthreads; it++) {
		nm2[i] += priv_nm2[it][i];
	    }
	}
    }

    ParallelDescriptor::ReduceRealSum(nm2.dataPtr(), n);

    for (int i=0; i<n; i++) {
	nm2[i] = std::sqrt(nm2[i]);
    }

    return nm2;
}
 
Real
MultiFab::norm1 (int comp, int ngrow, bool local) const
{
    Real nm1 = 0.e0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:nm1)
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        nm1 += get(mfi).norm(mfi.growntilebox(ngrow), 1, comp, 1);
    }

    if (!local)
	ParallelDescriptor::ReduceRealSum(nm1);

    return nm1;
}

Array<Real>
MultiFab::norm1 (const Array<int>& comps, int ngrow, bool local) const
{
    int n = comps.size();
    Array<Real> nm1(n, 0.e0);

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif
    PArray< Array<Real> > priv_nm1(nthreads, PArrayManage);
    for (int i=0; i<nthreads; i++) {
	priv_nm1.set(i, new Array<Real>(n, 0.0));
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
	int tid = omp_get_thread_num();
#else
	int tid = 0;
#endif
	for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
	{
            const Box& b = mfi.growntilebox(ngrow);
            for (int i=0; i<n; i++) {
                priv_nm1[tid][i] += get(mfi).norm(b, 1, comps[i], 1);
	    }
        }
#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
#endif
	for (int i=0; i<n; i++) {
	    for (int it=0; it<nthreads; it++) {
		nm1[i] += priv_nm1[it][i];
	    }
	}
    }

    if (!local)
	ParallelDescriptor::ReduceRealSum(nm1.dataPtr(), n);

    return nm1;
}

Real
MultiFab::sum (int comp, bool local) const
{
    Real sm = 0.e0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:sm)
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        sm += get(mfi).sum(mfi.tilebox(), comp, 1);
    }

    if (!local)
        ParallelDescriptor::ReduceRealSum(sm);

    return sm;
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
    BL_ASSERT(strt_comp + num_comp - 1 < n_comp && strt_comp + num_comp - 1 < mf.n_comp);
    BL_ASSERT(nghost <= n_grow && nghost <= mf.n_grow);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);

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
    BL_ASSERT(strt_comp + num_comp - 1 < n_comp && strt_comp + num_comp - 1 < mf.n_comp);
    BL_ASSERT(nghost <= n_grow && nghost <= mf.n_grow);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);

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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        get(mfi).plus(val,mfi.growntilebox(nghost),comp,num_comp);
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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        const Box& b = mfi.growntilebox(nghost) & region;

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
    BL_ASSERT(strt_comp + num_comp - 1 < n_comp && strt_comp + num_comp - 1 < mf.n_comp);
    BL_ASSERT(nghost <= n_grow && nghost <= mf.n_grow);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);

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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        get(mfi).mult(val, mfi.growntilebox(nghost), comp, num_comp);
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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        const Box& b = mfi.growntilebox(nghost) & region;

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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        get(mfi).invert(numerator, mfi.growntilebox(nghost), comp, num_comp);
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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        const Box& b = mfi.growntilebox(nghost) & region;

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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        get(mfi).negate(mfi.growntilebox(nghost), comp, num_comp);
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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
        const Box& b = mfi.growntilebox(nghost) & region;

        if (b.ok())
            get(mfi).negate(b,comp,num_comp);
    }
}

void
BoxLib::InterpAddBox (MultiFabCopyDescriptor& fabCopyDesc,
		      BoxList*                returnUnfilledBoxes,
		      Array<FillBoxId>&       returnedFillBoxIds,
		      const Box&              subbox,
		      MultiFabId              faid1,
		      MultiFabId              faid2,
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
BoxLib::InterpFillFab (MultiFabCopyDescriptor& fabCopyDesc,
		       const Array<FillBoxId>& fillBoxIds,
		       MultiFabId              faid1,
		       MultiFabId              faid2,
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
                       src_comp,
                       dest2,
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
                        bool local,
                        bool cross)
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

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	    std::vector< std::pair<int,Box> > isects;

	    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
	    {
		const int i = mfi.index();
		
		ba.intersections((*this)[mfi].box(),isects);
		
		for (int ii = 0, N = isects.size(); ii < N; ii++)
		{
		    const Box& bx  = isects[ii].second;
		    const int  iii = isects[ii].first;
		    
		    if (i != iii && DMap[iii] == MyProc)
		    {
			(*this)[mfi].copy((*this)[iii], bx, scomp, bx, scomp, ncomp);
		    }
		}
            }
        }
    }
    else
    {
        FabArray<FArrayBox>::FillBoundary(scomp,ncomp,cross);
    }
}

void
MultiFab::FillBoundary (bool local, bool cross)
{
    FillBoundary(0, n_comp, local, cross);
}

//
// Some useful typedefs.
//
typedef FabArrayBase::CopyComTag::CopyComTagsContainer CopyComTagsContainer;

typedef FabArrayBase::CopyComTag::MapOfCopyComTagContainers MapOfCopyComTagContainers;

void
MultiFab::SumBoundary (int scomp,
                       int ncomp)
{
    if ( n_grow <= 0 ) return;

    BL_PROFILE("MultiFab::SumBoundary()");
    //
    // We're going to attempt to reuse the information in the FillBoundary
    // cache.  The intersection info should be that same.  It's what we do
    // with it that's different.  Basically we have to send the m_RcvTags and
    // receive the m_SndTags, and invert the sense of fabIndex and srcIndex
    // in the CopyComTags.
    //
    MultiFab&                 mf       = *this;
    FabArrayBase::FBCacheIter cache_it = FabArrayBase::TheFB(false,mf);

    BL_ASSERT(cache_it != FabArrayBase::m_TheFBCache.end());

    const FabArrayBase::SI& TheSI = cache_it->second;

    if (ParallelDescriptor::NProcs() == 1)
    {
        //
        // There can only be local work to do.
        //
	int N_loc = (*TheSI.m_LocTags).size();
	// undafe to do OMP
	for (int i=0; i<N_loc; ++i)
        {
            const CopyComTag& tag = (*TheSI.m_LocTags)[i];
            mf[tag.srcIndex].plus(mf[tag.fabIndex],tag.box,tag.box,scomp,scomp,ncomp);
        }

        return;
    }

#ifdef BL_USE_MPI
    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    const int SeqNum = ParallelDescriptor::SeqNum();

    if (TheSI.m_LocTags->empty() && TheSI.m_RcvTags->empty() && TheSI.m_SndTags->empty())
        //
        // No work to do.
        //
        return;

    Array<MPI_Status>  stats;
    Array<int>         recv_from;
    Array<Real*>       recv_data;
    Array<MPI_Request> recv_reqs;
    //
    // Post rcvs. Allocate one chunk of space to hold'm all.
    //
    Real* the_recv_data = 0;

    FabArrayBase::PostRcvs(*TheSI.m_SndTags,*TheSI.m_SndVols,the_recv_data,recv_data,recv_from,recv_reqs,ncomp,SeqNum);

    //
    // Post send's
    //
    const int N_snds = TheSI.m_SndTags->size();

    Array<Real*>                       send_data;
    Array<int>                         send_N;
    Array<int>                         send_rank;
    Array<const CopyComTagsContainer*> send_cctc;

    send_data.reserve(N_snds);
    send_N   .reserve(N_snds);
    send_rank.reserve(N_snds);
    send_cctc.reserve(N_snds);

    for (MapOfCopyComTagContainers::const_iterator m_it = TheSI.m_RcvTags->begin(),
             m_End = TheSI.m_RcvTags->end();
         m_it != m_End;
         ++m_it)
    {
        std::map<int,int>::const_iterator vol_it = TheSI.m_RcvVols->find(m_it->first);

        BL_ASSERT(vol_it != TheSI.m_RcvVols->end());

        const int N = vol_it->second*ncomp;

        BL_ASSERT(N < std::numeric_limits<int>::max());

        Real* data = static_cast<Real*>(BoxLib::The_Arena()->alloc(N*sizeof(Real)));

	send_data.push_back(data);
	send_N   .push_back(N);
	send_rank.push_back(m_it->first);
	send_cctc.push_back(&(m_it->second));
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<N_snds; ++i)
    {
	Real*dptr = send_data[i];
	BL_ASSERT(dptr != 0);

	const CopyComTagsContainer& cctc = *send_cctc[i];

        for (CopyComTagsContainer::const_iterator it = cctc.begin();
             it != cctc.end(); ++it)
        {
            BL_ASSERT(distributionMap[it->fabIndex] == ParallelDescriptor::MyProc());
            const Box& bx = it->box;
            mf[it->fabIndex].copyToMem(bx,scomp,ncomp,dptr);
            const int Cnt = bx.numPts()*ncomp;
            dptr += Cnt;
        }
    }

    Array<MPI_Request> send_reqs;

    if (FabArrayBase::do_async_sends)
    {
	send_reqs.reserve(N_snds);
	for (int i=0; i<N_snds; ++i) {
            send_reqs.push_back(ParallelDescriptor::Asend
				(send_data[i],send_N[i],send_rank[i],SeqNum).req());
        }
    } else {
	for (int i=0; i<N_snds; ++i) {
            ParallelDescriptor::Send(send_data[i],send_N[i],send_rank[i],SeqNum);
            BoxLib::The_Arena()->free(send_data[i]);
        }
    }

    //
    // Do the local work.  Hope for a bit of communication/computation overlap.
    //
    int N_loc = (*TheSI.m_LocTags).size();
    // undafe to do OMP
    for (int i=0; i<N_loc; ++i)
    {
        const CopyComTag& tag = (*TheSI.m_LocTags)[i];

        BL_ASSERT(distributionMap[tag.fabIndex] == ParallelDescriptor::MyProc());
        BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

        mf[tag.srcIndex].plus(mf[tag.fabIndex],tag.box,tag.box,scomp,scomp,ncomp);
    }

    //
    //  wait and unpack
    //
    const int N_rcvs = TheSI.m_RcvTags->size();

    if (N_rcvs > 0)
    {
	Array<const CopyComTagsContainer*> recv_cctc;
	recv_cctc.reserve(N_rcvs);

        for (int k = 0; k < N_rcvs; k++)
        {
            MapOfCopyComTagContainers::const_iterator m_it = TheSI.m_SndTags->find(recv_from[k]);
            BL_ASSERT(m_it != TheSI.m_SndTags->end());

	    recv_cctc.push_back(&(m_it->second));
	}

	stats.resize(N_rcvs);
	BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, recv_reqs.dataPtr(), stats.dataPtr()) );

	// unsafe to do OMP
	{
	    FArrayBox fab;

	    for (int k = 0; k < N_rcvs; k++) 
	    {
		const Real* dptr = recv_data[k];
		BL_ASSERT(dptr != 0);
		
		const CopyComTagsContainer& cctc = *recv_cctc[k];
		
		for (CopyComTagsContainer::const_iterator it = cctc.begin();
		     it != cctc.end(); ++it)
		{
		    BL_ASSERT(distributionMap[it->srcIndex] == ParallelDescriptor::MyProc());
		    const Box& bx = it->box;
		    fab.resize(bx,ncomp);
		    const int Cnt = bx.numPts()*ncomp;
		    memcpy(fab.dataPtr(), dptr, Cnt*sizeof(Real));
		    mf[it->srcIndex].plus(fab,bx,bx,0,scomp,ncomp);
		    dptr += Cnt;
		}
	    }
        }
    }

    BoxLib::The_Arena()->free(the_recv_data);

    if (FabArrayBase::do_async_sends && !TheSI.m_RcvTags->empty())
        FabArrayBase::GrokAsyncSends(TheSI.m_RcvTags->size(),send_reqs,send_data,stats);

#endif /*BL_USE_MPI*/
}

void
MultiFab::SumBoundary ()
{
    SumBoundary(0, n_comp);
}


// Given a MultiFab in the compute MPI group, clone its data onto a MultiFab in
// the sidecar group.

// Note that the compute MultiFab will be a null pointer on the sidecar nodes,
// and the sidecar MultiFab will be null on the compute nodes. So be mindful of
// which processes will be executing which code when you access these pointers.
#ifdef BL_USE_MPI
void
MultiFab::SendMultiFabToSidecars (MultiFab *mf)
{
    // Broadcasts to intercommunicators have weird syntax. See below.
    // Whichever proc is actually broadcasting the data uses MPI_ROOT; all
    // other procs in that same communicator use MPI_PROC_NULL. On the
    // receiving end, the source rank is the rank of the broadcasting process
    // within its local communicator, *not* its rank in MPI_COMM_WORLD.
    // (Usually broadcasts will come from IOProcessor(), which is typically
    // rank 0.)
    const int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;

    if (ParallelDescriptor::Communicator() == ParallelDescriptor::CommunicatorComp())
    {
      // The in-transit procs also need the distribution map of the compute
      // group so they know who to MPI_Recv() from.
      Array<int> comp_procmap = mf->DistributionMap().ProcessorMap();
      int procmap_size = comp_procmap.size();
      ParallelDescriptor::Bcast(&procmap_size, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
      ParallelDescriptor::Bcast(&comp_procmap[0], procmap_size, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());

      // Even though a MultiFab is a distributed object, the layout of the
      // Boxes themselves is stored on every process, so broadcast them to
      // everyone on the in-transit side.
      const BoxArray& boxArray = mf->boxArray();
      for (int i = 0; i < boxArray.size(); ++i) {
        const Box& box = boxArray[i];
        const int *box_index_type = box.type().getVect();
        const int *smallEnd = box.smallEnd().getVect();
        const int *bigEnd = box.bigEnd().getVect();
        // getVect() requires a constant pointer, but MPI buffers require
        // non-constant pointers. Sorry this is awful.
        ParallelDescriptor::Bcast(const_cast<int*>(box_index_type), BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
        ParallelDescriptor::Bcast(const_cast<int*>(smallEnd)      , BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
        ParallelDescriptor::Bcast(const_cast<int*>(bigEnd)        , BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
      }

      int nComp = mf->nComp();
      ParallelDescriptor::Bcast(&nComp, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());

      int nGhost = mf->nGrow();
      ParallelDescriptor::Bcast(&nGhost, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());

      // Now that the sidecars have all the Box data, they will build the DM
      // and send it back to the compute processes, The compute DM has the same
      // size since they contain the same number of boxes.
      Array<int> intransit_procmap(mf->DistributionMap().ProcessorMap().size()+1);
      ParallelDescriptor::Bcast(&intransit_procmap[0], intransit_procmap.size(), 0, ParallelDescriptor::CommunicatorInter());

      // Both the compute and sidecar procs now have both DMs, so now we can
      // send the FAB data.
      for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
      {
          const int index = mfi.index();
          if (comp_procmap[index] == ParallelDescriptor::MyProc())
          {
            const Box& box = (*mf)[mfi].box();
            const long numPts = box.numPts();
            ParallelDescriptor::Send(&numPts, 1, intransit_procmap[index], 0, ParallelDescriptor::CommunicatorInter());
            const Real *dataPtr = (*mf)[mfi].dataPtr();
            ParallelDescriptor::Send(dataPtr, numPts*nComp, intransit_procmap[index], 1, ParallelDescriptor::CommunicatorInter());
          }
      }
    }
    else
    {
      // The sidecars also need the compute proc's DM.
      int procmap_size;
      ParallelDescriptor::Bcast(&procmap_size, 1, 0, ParallelDescriptor::CommunicatorInter());
      Array<int> comp_procmap(procmap_size);
      ParallelDescriptor::Bcast(&comp_procmap[0], procmap_size, 0, ParallelDescriptor::CommunicatorInter());
      DistributionMapping CompDM(comp_procmap);

      // Now build the list of Boxes.
      BoxList bl;

      for (unsigned int i = 0; i < procmap_size-1; ++i) {
        int box_index_type[BL_SPACEDIM];
        int smallEnd[BL_SPACEDIM];
        int bigEnd[BL_SPACEDIM];
        ParallelDescriptor::Bcast(box_index_type, BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());
        ParallelDescriptor::Bcast(smallEnd      , BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());
        ParallelDescriptor::Bcast(bigEnd        , BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());

        IntVect smallEnd_IV(smallEnd);
        IntVect bigEnd_IV(bigEnd);
        IntVect box_index_type_IV(box_index_type);
        Box box(smallEnd_IV, bigEnd_IV, box_index_type_IV);
        bl.push_back(box);
      }

      BoxArray ba(bl);

      // Get number of components in the FabArray.
      int nComp;
      ParallelDescriptor::Bcast(&nComp, 1, 0, ParallelDescriptor::CommunicatorInter());

      // Get number of ghost cells.
      int nGhost;
      ParallelDescriptor::Bcast(&nGhost, 1, 0, ParallelDescriptor::CommunicatorInter());

      // Now that the sidecars have all the Boxes, they can build their DM.
      DistributionMapping sidecar_DM;
      sidecar_DM.define(ba, ParallelDescriptor::NProcsSidecar());

      // The compute procs need the sidecars' DM so that we can match Send()s
      // and Recv()s for the FAB data.
      Array<int> intransit_procmap = sidecar_DM.ProcessorMap();
      ParallelDescriptor::Bcast(&intransit_procmap[0], procmap_size, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());

      // Now build the sidecar MultiFab with the new DM.
      mf->define(ba, nComp, nGhost, sidecar_DM, Fab_allocate);

      // Now we populate the MultiFab with data.
      for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
      {
          const int index = mfi.index();
        if (intransit_procmap[index] == ParallelDescriptor::MyProc())
        {
          long numPts;
          ParallelDescriptor::Recv(&numPts, 1, comp_procmap[index], 0, ParallelDescriptor::CommunicatorInter());
          Real FAB_data[numPts*nComp];
          Real *data_ptr = (*mf)[mfi].dataPtr();
          ParallelDescriptor::Recv(FAB_data, numPts*nComp, comp_procmap[index], 1, ParallelDescriptor::CommunicatorInter());
          std::memcpy(data_ptr, FAB_data, numPts*nComp*sizeof(Real));
        }
      }
    }
}
#endif


// If we want to use the multigrid solver from HPGMG then we must convert our
// MultiFabs to HPGMG's level data structures. This function essentially
// replaces the create_level() function in HPGMG.
#ifdef USEHPGMG
void
MultiFab::CreateHPGMGLevel (level_type* level,
                            const MultiFab& mf,
                            const int n_cell,
                            const int max_grid_size,
                            const int my_rank,
                            const int num_ranks,
                            const int domain_boundary_condition,
                            const int numVectors,
                            const double h0)
{
    int box;
    const int boxes_in_i = n_cell / max_grid_size;
    int TotalBoxes = boxes_in_i * boxes_in_i * boxes_in_i;

    // HPGMG requires perfect cubes for all boxes
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        if (!bx.isSquare()) {
             BoxLib::Error("All boxes must be square in HPGMG");
        }
    }

    // HPGMG also requires all boxes to be the same size, so we iterate over
    // all boxes and make sure they're the same.
    for (MFIter mfi1(mf); mfi1.isValid(); ++mfi1)
    {
        const Box& bx1 = mfi1.validbox();
        for (MFIter mfi2(mf); mfi2.isValid(); ++mfi2)
        {
            const Box& bx2 = mfi2.validbox();
            if (!(bx1.sameSize(bx2)))
            {
                BoxLib::Error("All boxes must be identical in HPGMG!");
            }
        }
    }

    // All the boxes have identical size and shape, so we just pick one of them
    // as a representative to fill in all the level data for HPGMG.
    MFIter mfi(mf);
    while (!mfi.isValid()) ++mfi;

    const Box& bx = mfi.validbox();
    const int box_dim = bx.length(0); /* Since we've already checked that all boxes are the same size, we can just use the size from one of them here. */

    if (TotalBoxes / num_ranks == 0)
      BoxLib::Error("Must have at least one box per MPI task when using HPGMG");

    if (ParallelDescriptor::IOProcessor())
    {
      std::cout << std::endl << "attempting to create a " << box_dim*boxes_in_i << "^3 level from " << TotalBoxes << " x " << box_dim << "^3 boxes distributed among " << num_ranks << " tasks..." << std::endl;
      if (domain_boundary_condition==BC_DIRICHLET)
      {
        std::cout << "boundary condition = BC_DIRICHLET" << std::endl;
      }
      else if (domain_boundary_condition==BC_PERIODIC)
      {
        std::cout << "boundary condition = BC_PERIODIC" << std::endl;
      }
      else
      {
        BoxLib::Error("Unknown boundary condition supplied");
      }
    }

    int omp_threads = 1;

#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp master
      {
        omp_threads = omp_get_num_threads ();
      }
    }
#endif

    int box_ghosts = stencil_get_radius();

    level->box_dim        = box_dim;
    level->box_ghosts     = box_ghosts;
    level->numVectors     = 0; // no vectors have been allocated yet
    level->vectors_base   = NULL; // pointer returned by bulk malloc
    level->vectors        = NULL; // pointers to individual vectors
    level->boxes_in.i     = boxes_in_i;
    level->boxes_in.j     = boxes_in_i;
    level->boxes_in.k     = boxes_in_i;
    level->dim.i          = box_dim*level->boxes_in.i;
    level->dim.j          = box_dim*level->boxes_in.j;
    level->dim.k          = box_dim*level->boxes_in.k;
    level->active         = 1;
    level->my_rank        = my_rank;
    level->num_ranks      = num_ranks;
    level->boundary_condition.type = domain_boundary_condition;
    level->must_subtract_mean = -1;
    level->num_threads      = omp_threads;
    level->my_blocks        = NULL;
    level->num_my_blocks    = 0;
    level->allocated_blocks = 0;
    level->tag              = log2(level->dim.i);
    level->h                = h0;
    level->fluxes           = NULL;

    // allocate 3D array of integers to hold the MPI rank of the corresponding box and initialize to -1 (unassigned)
    level->rank_of_box = (int*)malloc(level->boxes_in.i*level->boxes_in.j*level->boxes_in.k*sizeof(int));
    if(level->rank_of_box==NULL)
        BoxLib::Error("malloc of level->rank_of_box failed");
    for(box=0;box<level->boxes_in.i*level->boxes_in.j*level->boxes_in.k;box++){level->rank_of_box[box]=-1;}  // -1 denotes that there is no actual box assigned to this region


    // Now convert our rank distribution of boxes to HPGMG's rank_of_box array.
    // This is convoluted because HPGMG first assigns boxes to ranks, and then
    // lexicographically assigns the coordinates of each box. This
    // lexicographical ordering of box coordinates is *required* in order for
    // the MPI communication patterns in HPGMG to function correctly, via the
    // global_box_id variable. In other words, HPGMG anticipates the geometric
    // relationship between boxes based on their respective values of
    // global_box_id, and routes MPI traffic accordingly. However, in BoxLib
    // the box ranks and indices are not necessarily in this order, so we have
    // to "fake" the box ordering in HPGMG here (even though the coordinates
    // aren't actually assigned until we call create_vectors()) in order to
    // match the box ranks between BoxLib and HPGMG. This whole method is dumb
    // and deserves a better solution, but I don't know a better way to do it.

    int num_local_boxes = 0;
    int i,j,k;
    for(k=0;k<level->boxes_in.k;k++){
    for(j=0;j<level->boxes_in.j;j++){
    for(i=0;i<level->boxes_in.i;i++){
      int jStride = level->boxes_in.i;
      int kStride = level->boxes_in.i*level->boxes_in.j;
      int b=i + j*jStride + k*kStride;

      // These will be the coordinates of a box in HPGMG. These are also the
      // coordinates of a box already created in BoxLib. Now we iterate through
      // every rank's local boxes until we find the matching one, and assign
      // the rank of the HPGMG box to the same rank in BoxLib.

      const int low_i      = i*level->box_dim;
      const int low_j      = j*level->box_dim;
      const int low_k      = k*level->box_dim;

      bool found = false;
      for (MFIter mfi(mf); mfi.isValid(); ++mfi)
      {
        const Box &bx = mfi.validbox();
        const int *loVect = bx.loVect();

        // Found the matching box!
        if ((low_i == loVect[0]) &&
            (low_j == loVect[1]) &&
            (low_k == loVect[2]))
        {
            found = true;
            num_local_boxes++;
            break;
        }
      }
      if (found)
      {
        level->rank_of_box[b] = my_rank;
      }
    }}}

    // Now tell all the ranks what each other's box ranks are.
    const int tot_num_boxes = level->boxes_in.i * level->boxes_in.j * level->boxes_in.k;
    int all_box_ranks[tot_num_boxes];
    std::fill_n(all_box_ranks, tot_num_boxes, 1);
    MPI_Allreduce(level->rank_of_box, all_box_ranks, tot_num_boxes, MPI_INT, MPI_PROD, ParallelDescriptor::Communicator());
    for (unsigned int i = 0; i < tot_num_boxes; ++i)
    {
        level->rank_of_box[i] = std::abs(all_box_ranks[i]);
    }

    std::vector<int> box_ranks(level->rank_of_box, level->rank_of_box + tot_num_boxes);

    // calculate how many boxes I own...
    level->num_my_boxes=0;
    for(box=0;box<level->boxes_in.i*level->boxes_in.j*level->boxes_in.k;box++){if(level->rank_of_box[box]==level->my_rank)level->num_my_boxes++;}
    level->my_boxes = (box_type*)malloc(level->num_my_boxes*sizeof(box_type));
    if((level->num_my_boxes>0)&&(level->my_boxes==NULL))
        BoxLib::Error("malloc failed - create_level/level->my_boxes");

    // allocate flattened vector FP data and create pointers...
    if (ParallelDescriptor::IOProcessor())
        std::cout << "Allocating vectors... ";
    create_vectors (level, numVectors);
    if (ParallelDescriptor::IOProcessor())
        std::cout << "done." << std::endl;

    // Build and auxilarlly data structure that flattens boxes into blocks...
    for(box=0;box<level->num_my_boxes;box++){
      int blockcopy_i = BLOCKCOPY_TILE_I;
      int blockcopy_j = BLOCKCOPY_TILE_J;
      int blockcopy_k = BLOCKCOPY_TILE_K;

      append_block_to_list(&(level->my_blocks),&(level->allocated_blocks),&(level->num_my_blocks),
        /* dim.i         = */ level->my_boxes[box].dim,
        /* dim.j         = */ level->my_boxes[box].dim,
        /* dim.k         = */ level->my_boxes[box].dim,
        /* read.box      = */ box,
        /* read.ptr      = */ NULL,
        /* read.i        = */ 0,
        /* read.j        = */ 0,
        /* read.k        = */ 0,
        /* read.jStride  = */ level->my_boxes[box].jStride,
        /* read.kStride  = */ level->my_boxes[box].kStride,
        /* read.scale    = */ 1,
        /* write.box     = */ box,
        /* write.ptr     = */ NULL,
        /* write.i       = */ 0,
        /* write.j       = */ 0,
        /* write.k       = */ 0,
        /* write.jStride = */ level->my_boxes[box].jStride,
        /* write.kStride = */ level->my_boxes[box].kStride,
        /* write.scale   = */ 1,
        /* blockcopy_i   = */ blockcopy_i,
        /* blockcopy_j   = */ blockcopy_j,
        /* blockcopy_k   = */ blockcopy_k,
        /* subtype       = */ 0
      );
    }

    // build an assist structure for Gauss Seidel Red Black that would facilitate unrolling and SIMDization...
    level->RedBlack_base = NULL;
    level->RedBlack_FP = NULL;
    if(level->num_my_boxes){
      int i,j;
      int kStride = level->my_boxes[0].kStride;
      int jStride = level->my_boxes[0].jStride;
      level->RedBlack_base = (double*)malloc(2*kStride*sizeof(double)+256); // used for free()
      level->RedBlack_FP   = level->RedBlack_base; // aligned version
      // align first *non-ghost* zone element to a 64-Byte boundary...
      while( (uint64_t)(level->RedBlack_FP + level->box_ghosts*(1+level->box_jStride)) & 0x3f ){level->RedBlack_FP++;}
      // initialize RedBlack array...
      for(j=0-level->box_ghosts;j<level->box_dim+level->box_ghosts;j++){
      for(i=0-level->box_ghosts;i<level->box_dim+level->box_ghosts;i++){
        int ij = (i+level->box_ghosts) + (j+level->box_ghosts)*jStride;
        if((i^j^1)&0x1){
          level->RedBlack_FP[ij        ]=1.0;
          level->RedBlack_FP[ij+kStride]=0.0;
        }else{
          level->RedBlack_FP[ij        ]=0.0;
          level->RedBlack_FP[ij+kStride]=1.0;
        }
      }}
    }

    int shape;
    // create mini program for each stencil shape to perform a ghost zone exchange...
    for(shape=0;shape<STENCIL_MAX_SHAPES;shape++)build_exchange_ghosts(    level,shape);
    // create mini program for each stencil shape to perform a boundary condition...
    for(shape=0;shape<STENCIL_MAX_SHAPES;shape++)build_boundary_conditions(level,shape);


    // duplicate the parent communicator to be the communicator for each level
    #ifdef BL_USE_MPI
    if (ParallelDescriptor::IOProcessor())
        std::cout << "Duplicating MPI communicator... ";
    double time_start = MPI_Wtime();
    MPI_Comm_dup(ParallelDescriptor::Communicator(),&level->MPI_COMM_ALLREDUCE);
    double time_end = MPI_Wtime();
    double time_in_comm_dup = 0;
    double time_in_comm_dup_send = time_end-time_start;
    MPI_Allreduce(&time_in_comm_dup_send,&time_in_comm_dup,1,MPI_DOUBLE,MPI_MAX,ParallelDescriptor::Communicator());
    if (ParallelDescriptor::IOProcessor())
      std::cout << "done (" << time_in_comm_dup << " seconds)" << std::endl;
    #endif /* BL_USE_MPI */

    // report on potential load imbalance
    int BoxesPerProcess = level->num_my_boxes;
    #ifdef BL_USE_MPI
    int BoxesPerProcessSend = level->num_my_boxes;
    MPI_Allreduce(&BoxesPerProcessSend,&BoxesPerProcess,1,MPI_INT,MPI_MAX,ParallelDescriptor::Communicator());
    #endif /* BL_USE_MPI */
    if (ParallelDescriptor::IOProcessor())
      std::cout << "Calculating boxes per process... target=" << (double)TotalBoxes/(double)num_ranks << ", max=" << BoxesPerProcess << std::endl;
}


void MultiFab::SetupHPGMGCoefficients(const double a,
                                      const double b,
                                      const MultiFab& alpha,
                                      const MultiFab& beta_cc,
                                      level_type* level)
{

    // First set the alphas (cell-centered).
    bool found = false;
    for (MFIter mfi(alpha); mfi.isValid(); ++mfi) {

      const Box &bx = mfi.validbox();

      const int *loVect = bx.loVect();
      unsigned int box;
      for (box = 0; box < level->num_my_boxes; ++box)
      {
        if ((level->my_boxes[box].low.i == loVect[0]) &&
            (level->my_boxes[box].low.j == loVect[1]) &&
            (level->my_boxes[box].low.k == loVect[2]))
        {
          found = true;
          break;
        }
      }
      if (!found)
      {
        BoxLib::Error("Could not find matching boxes between HPGMG and BoxLib");
      }

      const Box &fabbox = mfi.fabbox();
      const double *alpha_data_ptr = alpha[mfi].dataPtr();
      int i,j,k;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int  ghosts = level->my_boxes[box].ghosts;
      const int   dim_i = level->my_boxes[box].dim;
      const int   dim_j = level->my_boxes[box].dim;
      const int   dim_k = level->my_boxes[box].dim;

      const int BL_jStride = fabbox.length(0);
      const int BL_kStride = fabbox.length(0) * fabbox.length(1);
      const int BoxLib_ghosts = alpha.nGrow();
      #ifdef _OPENMP
      #pragma omp parallel for private(k,j,i) collapse(3)
      #endif
      for(k=0;k<dim_k;k++){
      for(j=0;j<dim_j;j++){
      for(i=0;i<dim_i;i++){
        int ijk_HPGMG = (i+ghosts) + (j+ghosts)*jStride + (k+ghosts)*kStride;
        const int ijk_BoxLib = (i+BoxLib_ghosts) + (j+BoxLib_ghosts)*BL_jStride + (k+BoxLib_ghosts)*BL_kStride;
        level->my_boxes[box].vectors[VECTOR_ALPHA][ijk_HPGMG] = alpha_data_ptr[ijk_BoxLib];
      }}}
    }


    // Now convert the cell-centered beta to faces.
    found = false;
    for (MFIter mfi(beta_cc); mfi.isValid(); ++mfi) {

      const Box &bx = mfi.validbox();

      const int *loVect = bx.loVect();
      unsigned int box;
      for (box = 0; box < level->num_my_boxes; ++box)
      {
        if ((level->my_boxes[box].low.i == loVect[0]) &&
            (level->my_boxes[box].low.j == loVect[1]) &&
            (level->my_boxes[box].low.k == loVect[2]))
        {
          found = true;
          break;
        }
      }
      if (!found)
      {
        BoxLib::Error("Could not find matching boxes between HPGMG and BoxLib");
      }

      const Box &fabbox = mfi.fabbox();

      const double *beta_data_ptr = beta_cc[mfi].dataPtr();
      int i,j,k;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int  ghosts = level->my_boxes[box].ghosts;
      const int   dim_i = level->my_boxes[box].dim;
      const int   dim_j = level->my_boxes[box].dim;
      const int   dim_k = level->my_boxes[box].dim;
      const int BL_jStride = fabbox.length(0);
      const int BL_kStride = fabbox.length(0) * fabbox.length(1);
      const int BoxLib_ghosts = beta_cc.nGrow();

      #ifdef _OPENMP
      #pragma omp parallel for private(k,j,i) collapse(3)
      #endif
      for(k=0;k<=dim_k;k++){ // include high face
      for(j=0;j<=dim_j;j++){ // include high face
      for(i=0;i<=dim_i;i++){ // include high face
        //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        int ijk_HPGMG = (i+ghosts) + (j+ghosts)*jStride + (k+ghosts)*kStride;
        const int ijk_BoxLib   = (i  +BoxLib_ghosts) + (j  +BoxLib_ghosts)*BL_jStride + (k  +BoxLib_ghosts)*BL_kStride;
        const int im1jk_BoxLib = (i-1+BoxLib_ghosts) + (j  +BoxLib_ghosts)*BL_jStride + (k  +BoxLib_ghosts)*BL_kStride;
        const int ijm1k_BoxLib = (i  +BoxLib_ghosts) + (j-1+BoxLib_ghosts)*BL_jStride + (k  +BoxLib_ghosts)*BL_kStride;
        const int ijkm1_BoxLib = (i  +BoxLib_ghosts) + (j  +BoxLib_ghosts)*BL_jStride + (k-1+BoxLib_ghosts)*BL_kStride;
        level->my_boxes[box].vectors[VECTOR_BETA_I][ijk_HPGMG] = 0.5 * (beta_data_ptr[ijk_BoxLib] + beta_data_ptr[im1jk_BoxLib]);
        level->my_boxes[box].vectors[VECTOR_BETA_J][ijk_HPGMG] = 0.5 * (beta_data_ptr[ijk_BoxLib] + beta_data_ptr[ijm1k_BoxLib]);
        level->my_boxes[box].vectors[VECTOR_BETA_K][ijk_HPGMG] = 0.5 * (beta_data_ptr[ijk_BoxLib] + beta_data_ptr[ijkm1_BoxLib]);
      }}}
    }
}


void
MultiFab::ConvertToHPGMGLevel (const MultiFab& mf,
                               const int n_cell,
                               const int max_grid_size,
                               level_type* level,
                               const int component_id)
{
    bool found = false;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {

      const Box &bx = mfi.validbox();

      // The local box indices are ordered differently in HPGMG and BoxLib. So
      // as a simple (but SLOW) hack, we just find the boxes with matching
      // lower indices.
      // TODO: make this box matching less hacky
      const int *loVect = bx.loVect();
      unsigned int box;
      for (box = 0; box < level->num_my_boxes; ++box)
      {
        if ((level->my_boxes[box].low.i == loVect[0]) &&
            (level->my_boxes[box].low.j == loVect[1]) &&
            (level->my_boxes[box].low.k == loVect[2]))
        {
          found = true;
          break;
        }
      }

      if (!found)
      {
        BoxLib::Error("Could not find matching boxes between HPGMG and BoxLib");
      }

      const Box &fabbox = mfi.fabbox();
      const int BL_jStride = fabbox.length(0);
      const int BL_kStride = fabbox.length(0) * fabbox.length(1);

      const double *fab_data = mf[mfi].dataPtr();
      int i,j,k;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int  ghosts = level->my_boxes[box].ghosts;
      const int   dim_i = level->my_boxes[box].dim;
      const int   dim_j = level->my_boxes[box].dim;
      const int   dim_k = level->my_boxes[box].dim;
      const int BoxLib_ghosts = mf.nGrow();

      #ifdef _OPENMP
      #pragma omp parallel for private(k,j,i) collapse(3)
      #endif
      for(k=0;k<dim_k;k++){
      for(j=0;j<dim_j;j++){
      for(i=0;i<dim_i;i++){

    // The HPGMG strides are padded to align memory and encourage SIMD-ization,
    // so they are different than the BoxLib strides.

        const int ijk_HPGMG = (i+ghosts) + (j+ghosts)*jStride + (k+ghosts)*kStride;
        const int ijk_BoxLib = (i+BoxLib_ghosts) + (j+BoxLib_ghosts)*BL_jStride + (k+BoxLib_ghosts)*BL_kStride;

        level->my_boxes[box].vectors[component_id][ijk_HPGMG] = fab_data[ijk_BoxLib];

      }}}

    }
}

void MultiFab::ConvertFromHPGMGLevel(MultiFab& mf,
                                     const level_type* level,
                                     const int component_id)
{
  for (MFIter mfi(mf); mfi.isValid(); ++mfi)
  {
      const Box &bx = mfi.validbox();
      double *fab_data = mf[mfi].dataPtr();

      // First find the HPGMG box corresponding to this BoxLib box.
      const int *loVect = bx.loVect();
      int box;
      for (box = 0; box < level->num_my_boxes; ++box)
      {
        if ((level->my_boxes[box].low.i == loVect[0]) &&
            (level->my_boxes[box].low.j == loVect[1]) &&
            (level->my_boxes[box].low.k == loVect[2]))
          break;
      }

      const Box &fabbox = mfi.fabbox();

      // Found the matching boxes, now fill the data.
      const int dim_i = level->my_boxes[box].dim;
      const int dim_j = level->my_boxes[box].dim;
      const int dim_k = level->my_boxes[box].dim;
      const int ghosts = level->my_boxes[box].ghosts;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int BoxLib_ghosts = mf.nGrow();

      int i, j, k;
      #ifdef _OPENMP
      #pragma omp parallel for private(k,j,i) collapse(3)
      #endif
      for(k=0;k<dim_k;k++){
      for(j=0;j<dim_j;j++){
      for(i=0;i<dim_i;i++){

        const int ijk_HPGMG = (i+ghosts) + (j+ghosts)*jStride + (k+ghosts)*kStride;

        // WARNING: this indexing stride works for FABs *ONLY* if we have ONE
        // component in the FAB. If we have more than one we have to stride
        // over the components in the outermost loop (outside of k).
        const int BL_jStride = fabbox.length(0);
        const int BL_kStride = fabbox.length(0) * fabbox.length(1);
        const int ijk_BoxLib = (i+BoxLib_ghosts) + (j+BoxLib_ghosts)*BL_jStride + (k+BoxLib_ghosts)*BL_kStride;

        fab_data[ijk_BoxLib] = level->my_boxes[box].vectors[VECTOR_U][ijk_HPGMG];
      }}}
  }
}
#endif /* USEHPGMG */
