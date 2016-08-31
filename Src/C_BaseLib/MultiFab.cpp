
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

#ifdef BL_MEM_PROFILING
#include <MemProfiler.H>
#endif

namespace
{
    bool initialized = false;
    int num_multifabs     = 0;
    int num_multifabs_hwm = 0;
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
MultiFab::Saxpy (MultiFab&       dst,
		 Real            a, 
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
            dst[mfi].saxpy(a, src[mfi], bx, bx, srccomp, dstcomp, numcomp);
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

    FArrayBox foo(); // This will call FArrayBox::Initialize().

    BoxLib::ExecOnFinalize(MultiFab::Finalize);

#ifdef BL_MEM_PROFILING
    MemProfiler::add("MultiFab", std::function<MemProfiler::NBuildsInfo()>
		     ([] () -> MemProfiler::NBuildsInfo {
			 return {num_multifabs, num_multifabs_hwm};
		     }));
#endif

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
#ifdef BL_MEM_PROFILING
    ++num_multifabs;
    num_multifabs_hwm = std::max(num_multifabs_hwm, num_multifabs);
#endif
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
    if (SharedMemory() && alloc == Fab_allocate) initVal();  // else already done in FArrayBox
#ifdef BL_MEM_PROFILING
    ++num_multifabs;
    num_multifabs_hwm = std::max(num_multifabs_hwm, num_multifabs);
#endif
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
    if (SharedMemory() && alloc == Fab_allocate) initVal();  // else already done in FArrayBox
#ifdef BL_MEM_PROFILING
    ++num_multifabs;
    num_multifabs_hwm = std::max(num_multifabs_hwm, num_multifabs);
#endif
}

MultiFab::MultiFab (const BoxArray& bxs,
                    int             ncomp,
                    int             ngrow,
		    ParallelDescriptor::Color color)
    :
    FabArray<FArrayBox>(bxs,ncomp,ngrow,color)
{
    Initialize();
    if (SharedMemory()) initVal();  // else already done in FArrayBox
#ifdef BL_MEM_PROFILING
    ++num_multifabs;
    num_multifabs_hwm = std::max(num_multifabs_hwm, num_multifabs);
#endif
}

MultiFab::~MultiFab()
{
#ifdef BL_MEM_PROFILING
    --num_multifabs;
#endif    
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
		  const IntVect&  nodal,
		  ParallelDescriptor::Color color)
{
    this->FabArray<FArrayBox>::define(bxs,nvar,ngrow,alloc,nodal,color);
    if (SharedMemory() && alloc == Fab_allocate) initVal();  // else already done in FArrayBox
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
    if (SharedMemory() && alloc == Fab_allocate) initVal();  // else already done in FArrayBox
}

void
MultiFab::initVal ()
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
	(*this)[mfi].initVal();
    }
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
	ParallelDescriptor::ReduceBoolOr(r,this->color());

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
	ParallelDescriptor::ReduceBoolOr(r,this->color());

    return r;
}

bool 
MultiFab::contains_inf (bool local) const
{
    return contains_inf(0,nComp(),nGrow(),local);
}

bool 
MultiFab::is_nodal () const
{
    return boxArray().ixType().nodeCentered();
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
	ParallelDescriptor::ReduceRealMin(mn,this->color());

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
	ParallelDescriptor::ReduceRealMin(mn,this->color());

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
	ParallelDescriptor::ReduceRealMax(mx,this->color());

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
	ParallelDescriptor::ReduceRealMax(mx,this->color());

    return mx;
}

IntVect
MultiFab::minIndex (int comp,
                    int nghost) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(this->color() == ParallelDescriptor::DefaultColor());

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
    BL_ASSERT(this->color() == ParallelDescriptor::DefaultColor());

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
MultiFab::norm0 (int comp, const BoxArray& ba, int nghost, bool local) const
{
    Real nm0 = -std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel reduction(max:nm0)
#endif
    {
	std::vector< std::pair<int,Box> > isects;

	for (MFIter mfi(*this); mfi.isValid(); ++mfi)
	{
	    ba.intersections(BoxLib::grow(mfi.validbox(),nghost),isects);

	    for (int i = 0, N = isects.size(); i < N; i++)
	    {
		nm0 = std::max(nm0, get(mfi).norm(isects[i].second, 0, comp, 1));
	    }
	}
    }
 
    if (!local)
	ParallelDescriptor::ReduceRealMax(nm0,this->color());
 
    return nm0;
}

Real
MultiFab::norm0 (int comp, int nghost, bool local) const
{
    Real nm0 = -std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel reduction(max:nm0)
#endif
    for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
    {
	nm0 = std::max(nm0, get(mfi).norm(mfi.growntilebox(nghost), 0, comp, 1));
    }

    if (!local)
	ParallelDescriptor::ReduceRealMax(nm0,this->color());

    return nm0;
}

Array<Real>
MultiFab::norm0 (const Array<int>& comps, int nghost, bool local) const
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
	        priv_nm0[tid][i] = std::max(priv_nm0[tid][i], 
					    get(mfi).norm(mfi.growntilebox(nghost), 0, comps[i], 1));
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
	ParallelDescriptor::ReduceRealMax(nm0.dataPtr(), n, this->color());

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

    ParallelDescriptor::ReduceRealSum(nm2,this->color());

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

    ParallelDescriptor::ReduceRealSum(nm2.dataPtr(), n, this->color());

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
	ParallelDescriptor::ReduceRealSum(nm1,this->color());

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
	ParallelDescriptor::ReduceRealSum(nm1.dataPtr(), n, this->color());

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
        ParallelDescriptor::ReduceRealSum(sm, this->color());

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
	dest1.setVal(std::numeric_limits<Real>::quiet_NaN());
        FArrayBox dest2(dest.box(), dest.nComp());
	dest2.setVal(std::numeric_limits<Real>::quiet_NaN());
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

//
// Some useful typedefs.
//
typedef FabArrayBase::CopyComTag::CopyComTagsContainer CopyComTagsContainer;

typedef FabArrayBase::CopyComTag::MapOfCopyComTagContainers MapOfCopyComTagContainers;

void
MultiFab::SumBoundary (int scomp, int ncomp, const Periodicity& period)
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
    MultiFab& mf = *this;

    const FabArrayBase::FB& TheFB = mf.getFB(period);

    const CopyComTagsContainer&      LocTags = *(TheFB.m_LocTags);
    const MapOfCopyComTagContainers& SndTags = *(TheFB.m_RcvTags);
    const MapOfCopyComTagContainers& RcvTags = *(TheFB.m_SndTags);
    const std::map<int,int>&         SndVols = *(TheFB.m_RcvVols);
    const std::map<int,int>&         RcvVols = *(TheFB.m_SndVols);

    if (ParallelDescriptor::NProcs() == 1)
    {
        //
        // There can only be local work to do.
        //
	int N_loc = LocTags.size();
	// undafe to do OMP
	for (int i=0; i<N_loc; ++i)
        {
            const CopyComTag& tag = LocTags[i];
            mf[tag.srcIndex].plus(mf[tag.dstIndex],tag.dbox,tag.sbox,scomp,scomp,ncomp);
        }

        return;
    }

#ifdef BL_USE_MPI
    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int SeqNum;
    {
	ParallelDescriptor::Color mycolor = this->color();
	if (mycolor == ParallelDescriptor::DefaultColor()) {
	    SeqNum = ParallelDescriptor::SeqNum();
	} else if (mycolor == ParallelDescriptor::SubCommColor()) {
	    SeqNum = ParallelDescriptor::SubSeqNum();
	}
	// else I don't have any data and my SubSeqNum() should not be called.
    }

    if (LocTags.empty() && RcvTags.empty() && SndTags.empty())
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

    FabArrayBase::PostRcvs(RcvVols,the_recv_data,recv_data,recv_from,recv_reqs,ncomp,SeqNum);

    //
    // Post send's
    //
    const int N_snds = SndTags.size();

    Array<Real*>                       send_data;
    Array<int>                         send_N;
    Array<int>                         send_rank;
    Array<const CopyComTagsContainer*> send_cctc;

    send_data.reserve(N_snds);
    send_N   .reserve(N_snds);
    send_rank.reserve(N_snds);
    send_cctc.reserve(N_snds);

    for (MapOfCopyComTagContainers::const_iterator m_it = SndTags.begin(),
             m_End = SndTags.end();
         m_it != m_End;
         ++m_it)
    {
        std::map<int,int>::const_iterator vol_it = SndVols.find(m_it->first);

        BL_ASSERT(vol_it != SndVols.end());

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
            BL_ASSERT(distributionMap[it->dstIndex] == ParallelDescriptor::MyProc());
            const Box& bx = it->dbox;
            mf[it->dstIndex].copyToMem(bx,scomp,ncomp,dptr);
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
    int N_loc = LocTags.size();
    // undafe to do OMP
    for (int i=0; i<N_loc; ++i)
    {
        const CopyComTag& tag = LocTags[i];

        BL_ASSERT(distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc());
        BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

        mf[tag.srcIndex].plus(mf[tag.dstIndex],tag.dbox,tag.sbox,scomp,scomp,ncomp);
    }

    //
    //  wait and unpack
    //
    const int N_rcvs = RcvTags.size();

    if (N_rcvs > 0)
    {
	Array<const CopyComTagsContainer*> recv_cctc;
	recv_cctc.reserve(N_rcvs);

        for (int k = 0; k < N_rcvs; k++)
        {
            MapOfCopyComTagContainers::const_iterator m_it = RcvTags.find(recv_from[k]);
            BL_ASSERT(m_it != RcvTags.end());

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
		    const Box& bx = it->sbox;
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

    if (FabArrayBase::do_async_sends && !SndTags.empty())
        FabArrayBase::WaitForAsyncSends(N_snds,send_reqs,send_data,stats);

#endif /*BL_USE_MPI*/
}

void
MultiFab::SumBoundary (const Periodicity& period)
{
    SumBoundary(0, n_comp, period);
}


void
MultiFab::AddProcsToComp (int ioProcNumSCS, int ioProcNumAll,
                          int scsMyId, MPI_Comm scsComm)
{
  // ---- bools
  int bInit(initialized);
  if(scsMyId != ioProcNumSCS) {
    initialized   = bInit;
  }
  FabArray::AddProcsToComp(ioProcNumSCS, ioProcNumAll, scsMyId, scsComm);
}

