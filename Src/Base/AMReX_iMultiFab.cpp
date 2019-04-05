
#include <algorithm>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <map>
#include <limits>
#include <climits>

#include <AMReX_BLassert.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>

namespace amrex {

namespace
{
    bool initialized = false;
}

void
iMultiFab::Add (iMultiFab&       dst,
	       const iMultiFab& src,
	       int             srccomp,
	       int             dstcomp,
	       int             numcomp,
	       int             nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrow() >= nghost && src.nGrow() >= nghost);

    amrex::Add(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
iMultiFab::Copy (iMultiFab&       dst,
                const iMultiFab& src,
                int             srccomp,
                int             dstcomp,
                int             numcomp,
                int             nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrow() >= nghost && src.nGrow() >= nghost);

    amrex::Copy(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
iMultiFab::Subtract (iMultiFab&       dst,
		    const iMultiFab& src,
		    int             srccomp,
		    int             dstcomp,
		    int             numcomp,
		    int             nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrow() >= nghost && src.nGrow() >= nghost);

    amrex::Subtract(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
iMultiFab::Multiply (iMultiFab&       dst,
		    const iMultiFab& src,
		    int             srccomp,
		    int             dstcomp,
		    int             numcomp,
		    int             nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrow() >= nghost && src.nGrow() >= nghost);

    amrex::Multiply(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
iMultiFab::Divide (iMultiFab&       dst,
		  const iMultiFab& src,
		  int             srccomp,
		  int             dstcomp,
		  int             numcomp,
		  int             nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrow() >= nghost && src.nGrow() >= nghost);

    amrex::Divide(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
iMultiFab::plus (int val,
                 int  nghost)
{
    plus(val,0,n_comp,nghost);
}

void
iMultiFab::plus (int       val,
                 const Box& region,
                 int        nghost)
{
    plus(val,region,0,n_comp,nghost);
}

void
iMultiFab::mult (int val,
                 int  nghost)
{
    mult(val,0,n_comp,nghost);
}

void
iMultiFab::mult (int       val,
                 const Box& region,
                 int        nghost)
{
    mult(val,region,0,n_comp,nghost);
}

void
iMultiFab::negate (int nghost)
{
    negate(0,n_comp,nghost);
}

void
iMultiFab::negate (const Box& region,
                  int        nghost)
{
    negate(region,0,n_comp,nghost);
}

void
iMultiFab::Initialize ()
{
    if (initialized) return;

    amrex::ExecOnFinalize(iMultiFab::Finalize);

    initialized = true;
}

void
iMultiFab::Finalize ()
{
    initialized = false;
}

iMultiFab::iMultiFab () {}

iMultiFab::iMultiFab (const BoxArray&            bxs,
                      const DistributionMapping& dm,
                      int                        ncomp,
                      int                        ngrow,
		      const MFInfo&              info,
                      const FabFactory<IArrayBox>& factory)
    : iMultiFab(bxs,dm,ncomp,IntVect(ngrow),info,factory)
{
}

iMultiFab::iMultiFab (const BoxArray&            bxs,
                      const DistributionMapping& dm,
                      int                        ncomp,
                      const IntVect&             ngrow,
                      const MFInfo&              info,
                      const FabFactory<IArrayBox>& factory)
    :
    FabArray<IArrayBox>(bxs,dm,ncomp,ngrow,info,factory)
{
}

iMultiFab::iMultiFab (const iMultiFab& rhs, MakeType maketype, int scomp, int ncomp)
    :
    FabArray<IArrayBox>(rhs, maketype, scomp, ncomp)
{
}

void
iMultiFab::operator= (int r)
{
    setVal(r);
}

void
iMultiFab::define (const BoxArray&            bxs,
		   const DistributionMapping& dm,
		   int                        nvar,
		   const IntVect&             ngrow,
		   const MFInfo&              info,
                   const FabFactory<IArrayBox>& factory)
{
    this->FabArray<IArrayBox>::define(bxs,dm,nvar,ngrow,info, factory);
}

void
iMultiFab::define (const BoxArray&            bxs,
		   const DistributionMapping& dm,
		   int                        nvar,
		   int                        ngrow,
		   const MFInfo&              info,
                   const FabFactory<IArrayBox>& factory)
{
    this->FabArray<IArrayBox>::define(bxs,dm,nvar,ngrow,info, factory);
}

int
iMultiFab::min (int comp,
		int nghost,
		bool local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    int mn = amrex::ReduceMin(*this, nghost,
    [=] AMREX_GPU_HOST_DEVICE (Box const& bx, IArrayBox const& fab) -> int
    {
        return fab.min(bx,comp);
    });                          

    if (!local)
	ParallelDescriptor::ReduceIntMin(mn);

    return mn;
}

int
iMultiFab::min (const Box& region,
                int        comp,
                int        nghost,
		bool       local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    int mn = amrex::ReduceMin(*this, nghost,
    [=] AMREX_GPU_HOST_DEVICE (Box const& bx, IArrayBox const& fab) -> int
    {
        const Box& b = bx & region;
        if (b.ok()) {
            return fab.min(b,comp);
        } else {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
            return std::numeric_limits<int>::max();
#else
            return INT_MAX;
#endif
        }
    });

    if (!local)
	ParallelDescriptor::ReduceIntMin(mn);

    return mn;
}

int
iMultiFab::max (int comp,
		int nghost,
		bool local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    int mx = amrex::ReduceMax(*this, nghost,
    [=] AMREX_GPU_HOST_DEVICE (Box const& bx, IArrayBox const& fab) -> int
    {
        return fab.max(bx,comp);
    });                          

    if (!local)
	ParallelDescriptor::ReduceIntMax(mx);

    return mx;
}

int
iMultiFab::max (const Box& region,
		int        comp,
		int        nghost,
		bool       local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    int mx = amrex::ReduceMax(*this, nghost,
    [=] AMREX_GPU_HOST_DEVICE (Box const& bx, IArrayBox const& fab) -> int
    {
        const Box& b = bx & region;
        if (b.ok()) {
            return fab.max(b,comp);
        } else {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
            return std::numeric_limits<int>::lowest();
#else
            return INT_MIN;
#endif
        }
    });

    if (!local)
	ParallelDescriptor::ReduceIntMax(mx);

    return mx;
}

long
iMultiFab::sum (int comp, int nghost, bool local) const
{
    AMREX_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    iMultiFab imf(*this, amrex::make_alias, comp, 1);
    FabArray<BaseFab<long> > lmf = ToLongMultiFab(imf);

    long sm = amrex::ReduceSum(lmf, nghost,
    [=] AMREX_GPU_HOST_DEVICE (Box const& bx, BaseFab<long> const& fab) -> long
    {
        return fab.sum(bx,0);
    });

    if (!local) ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());

    return sm;
}

IntVect
iMultiFab::minIndex (int comp,
                    int nghost) const
{
    // TODO GPU

    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    IntVect loc;

    int mn = std::numeric_limits<int>::max();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	IntVect priv_loc;
	int priv_mn = std::numeric_limits<int>::max();
	
	for (MFIter mfi(*this); mfi.isValid(); ++mfi)
	{
	    const Box& bx = amrex::grow(mfi.validbox(),nghost);
	    const int  lmn = get(mfi).min(bx,comp);
	    
	    if (lmn < priv_mn)
	    {
		priv_mn  = lmn;
		priv_loc = get(mfi).minIndex(bx,comp);
	    }
	}

#ifdef _OPENMP
#pragma omp critical (imultifab_minindex)
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
        Vector<int> mns(1);
        Vector<int>  locs(1);

        if (ParallelDescriptor::IOProcessor())
        {
            mns.resize(NProcs);
            locs.resize(NProcs*AMREX_SPACEDIM);
        }

        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::Gather(&mn, 1, mns.dataPtr(), 1, IOProc);

        BL_ASSERT(sizeof(IntVect) == sizeof(int)*AMREX_SPACEDIM);

        ParallelDescriptor::Gather(loc.getVect(), AMREX_SPACEDIM, locs.dataPtr(), AMREX_SPACEDIM, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            mn  = mns[0];
            loc = IntVect(AMREX_D_DECL(locs[0],locs[1],locs[2]));

            for (int i = 1; i < NProcs; i++)
            {
                if (mns[i] < mn)
                {
                    mn = mns[i];

                    const int j = AMREX_SPACEDIM * i;

                    loc = IntVect(AMREX_D_DECL(locs[j+0],locs[j+1],locs[j+2]));
                }
            }
        }

        ParallelDescriptor::Bcast(const_cast<int*>(loc.getVect()), AMREX_SPACEDIM, IOProc);
    }

    return loc;
}

IntVect
iMultiFab::maxIndex (int comp,
                    int nghost) const
{
    // TODO

    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    IntVect loc;

    int mx = -std::numeric_limits<int>::max();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	IntVect priv_loc;
	int priv_mx = -std::numeric_limits<int>::max();

	for (MFIter mfi(*this); mfi.isValid(); ++mfi)
	{
	    const Box& bx = amrex::grow(mfi.validbox(),nghost);
	    const int  lmx = get(mfi).max(bx,comp);
	    
	    if (lmx > priv_mx)
	    {
		priv_mx  = lmx;
		priv_loc = get(mfi).maxIndex(bx,comp);
	    }
	}

#ifdef _OPENMP
#pragma omp critical (imultifab_maxindex)
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
        Vector<int> mxs(1);
        Vector<int>  locs(1);

        if (ParallelDescriptor::IOProcessor())
        {
            mxs.resize(NProcs);
            locs.resize(NProcs*AMREX_SPACEDIM);
        }

        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::Gather(&mx, 1, mxs.dataPtr(), 1, IOProc);

        BL_ASSERT(sizeof(IntVect) == sizeof(int)*AMREX_SPACEDIM);

        ParallelDescriptor::Gather(loc.getVect(), AMREX_SPACEDIM, locs.dataPtr(), AMREX_SPACEDIM, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            mx  = mxs[0];
            loc = IntVect(AMREX_D_DECL(locs[0],locs[1],locs[2]));

            for (int i = 1; i < NProcs; i++)
            {
                if (mxs[i] > mx)
                {
                    mx = mxs[i];

                    const int j = AMREX_SPACEDIM * i;

                    loc = IntVect(AMREX_D_DECL(locs[j+0],locs[j+1],locs[j+2]));
                }
            }
        }

        ParallelDescriptor::Bcast(const_cast<int*>(loc.getVect()), AMREX_SPACEDIM, IOProc);
    }

    return loc;
}

void
iMultiFab::minus (const iMultiFab& mf,
                 int             strt_comp,
                 int             num_comp,
                 int             nghost)
{
    iMultiFab::Subtract(*this, mf, strt_comp, strt_comp, num_comp, nghost);
}

void
iMultiFab::divide (const iMultiFab& mf,
		  int             strt_comp,
		  int             num_comp,
		  int             nghost)
{
    iMultiFab::Divide(*this, mf, strt_comp, strt_comp, num_comp, nghost);
}

void
iMultiFab::plus (int val,
                 int  comp,
                 int  num_comp,
                 int  nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    FabArray<IArrayBox>::plus(val,comp,num_comp,nghost);
}

void
iMultiFab::plus (int       val,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    FabArray<IArrayBox>::plus(val,region,comp,num_comp,nghost);
}

void
iMultiFab::plus (const iMultiFab& mf,
                int             strt_comp,
                int             num_comp,
                int             nghost)
{
    BL_ASSERT(boxarray == mf.boxarray);
    BL_ASSERT(strt_comp >= 0);
    BL_ASSERT(num_comp > 0);
    BL_ASSERT(strt_comp + num_comp - 1 < n_comp && strt_comp + num_comp - 1 < mf.n_comp);
    BL_ASSERT(nghost <= n_grow.min() && nghost <= mf.n_grow.min());

    amrex::Add(*this, mf, strt_comp, strt_comp, num_comp, nghost);
}

void
iMultiFab::mult (int val,
                 int  comp,
                 int  num_comp,
                 int  nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    FabArray<IArrayBox>::mult(val,comp,num_comp,nghost);
}

void
iMultiFab::mult (int       val,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    FabArray<IArrayBox>::mult(val,region,comp,num_comp,nghost);
}

void
iMultiFab::negate (int comp,
                  int num_comp,
                  int nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);

    FabArray<IArrayBox>::mult(-1,comp,num_comp,nghost);
}

void
iMultiFab::negate (const Box& region,
                  int        comp,
                  int        num_comp,
                  int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);

    FabArray<IArrayBox>::mult(-1,region,comp,num_comp,nghost);
}

std::unique_ptr<iMultiFab>
OwnerMask (FabArrayBase const& mf, const Periodicity& period)
{
    //TODO GPU????
    BL_PROFILE("OwnerMask()");

    const BoxArray& ba = mf.boxArray();
    const DistributionMapping& dm = mf.DistributionMap();

    const int owner = 1;
    const int nonowner = 0;

    std::unique_ptr<iMultiFab> p{new iMultiFab(ba,dm,1,0, MFInfo(),
                                               DefaultFabFactory<IArrayBox>())};
    p->setVal(owner);

    const std::vector<IntVect>& pshifts = period.shiftIntVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector< std::pair<int,Box> > isects;
        
        for (MFIter mfi(*p); mfi.isValid(); ++mfi)
        {
            IArrayBox& fab = (*p)[mfi];
            const Box& bx = fab.box();
            const int i = mfi.index();
            for (const auto& iv : pshifts)
            {
                ba.intersections(bx+iv, isects);                    
                for (const auto& is : isects)
                {
                    const int oi = is.first;
                    const Box& obx = is.second;
                    if ((oi < i) || (oi == i && iv < IntVect::TheZeroVector())) {
                        fab.setVal(nonowner, obx-iv, 0, 1);
                    }
                }
            }
        }
    }

    return p;
}

}
