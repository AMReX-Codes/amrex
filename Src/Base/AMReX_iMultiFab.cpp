
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

    long sm = 0;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<long> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(*this); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            const auto& arr = this->array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return { static_cast<long>(arr(i,j,k,comp)) };
            });
        }

        ReduceTuple hv = reduce_data.value();
        sm = amrex::get<0>(hv);
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            sm += (*this)[mfi].sum(bx,comp);
        }
    }

    if (!local) ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());

    return sm;
}

namespace {

static IntVect
indexFromValue (iMultiFab const& mf, int comp, int nghost, int value, MPI_Op mmloc)
{
    IntVect loc;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        int tmp[1+AMREX_SPACEDIM] = {0};
        amrex::Gpu::AsyncArray<int> aa(tmp, 1+AMREX_SPACEDIM);
        int* p = aa.data();
        // This is a device ptr to 1+AMREX_SPACEDIM int zeros.
        // The first is used as an atomic bool and the others for intvect.
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& bx = amrex::grow(mfi.validbox(), nghost);
            const Array4<int const> arr = mf.array(mfi);
            AMREX_LAUNCH_DEVICE_LAMBDA(bx, tbx,
            {
                int* flag = p;
                if (*flag == 0) {
                    const IArrayBox fab(arr);
                    IntVect t_loc = fab.indexFromValue(value, tbx, comp);
                    if (tbx.contains(t_loc)) {
                        if (Gpu::Atomic::Exch(flag,1) == 0) {
                            AMREX_D_TERM(p[1] = t_loc[0];,
                                         p[2] = t_loc[1];,
                                         p[3] = t_loc[2];);
                        }
                    }
                }
            });
        }
        aa.copyToHost(tmp, 1+AMREX_SPACEDIM);
        AMREX_D_TERM(loc[0] = tmp[1];,
                     loc[1] = tmp[2];,
                     loc[2] = tmp[3];);
    }
    else
#endif
    {
        bool f = false;
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            IntVect priv_loc = IntVect::TheMinVector();
            for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox(nghost);
                const IArrayBox& fab = mf[mfi];
                IntVect t_loc = fab.indexFromValue(value, bx, comp);
                if (bx.contains(t_loc)) {
                    priv_loc = t_loc;
                };
            }

            if (priv_loc.allGT(IntVect::TheMinVector())) {
                bool old;
// we should be able to test on _OPENMP < 201107 for capture (version 3.1)
// but we must work around a bug in gcc < 4.9
#if defined(_OPENMP) && _OPENMP < 201307
#pragma omp critical (amrex_indexfromvalue)
#elif defined(_OPENMP)
#pragma omp atomic capture
#endif
                {
                    old = f;
                    f = true;
                }

                if (old == false) loc = priv_loc;
            }
        }
    }

#ifdef BL_USE_MPI
    const int NProcs = ParallelContext::NProcsSub();
    if (NProcs > 1)
    {
        struct {
            int mm;
            int rank;
        } in, out;
        in.mm = value;
        in.rank = ParallelContext::MyProcSub();
        MPI_Datatype datatype = MPI_2INT;
        MPI_Comm comm = ParallelContext::CommunicatorSub();
        MPI_Allreduce(&in,  &out, 1, datatype, mmloc, comm);
        MPI_Bcast(&(loc[0]), AMREX_SPACEDIM, MPI_INT, out.rank, comm);
    }
#endif

    return loc;
}

}

IntVect
iMultiFab::minIndex (int comp, int nghost) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    int mn = this->min(comp, nghost, true);
    return indexFromValue(*this, comp, nghost, mn, MPI_MINLOC);
}

IntVect
iMultiFab::maxIndex (int comp, int nghost) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    int mx = this->max(comp, nghost, true);
    return indexFromValue(*this, comp, nghost, mx, MPI_MAXLOC);
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
    BL_PROFILE("OwnerMask()");

    const BoxArray& ba = mf.boxArray();
    const DistributionMapping& dm = mf.DistributionMap();

    const int owner = 1;
    const int nonowner = 0;

    std::unique_ptr<iMultiFab> p{new iMultiFab(ba,dm,1,0, MFInfo(),
                                               DefaultFabFactory<IArrayBox>())};
    const std::vector<IntVect>& pshifts = period.shiftIntVect();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        std::vector< std::pair<int,Box> > isects;
        
        for (MFIter mfi(*p); mfi.isValid(); ++mfi)
        {
            const Box& bx = (*p)[mfi].box();
            auto arr = p->array(mfi);
            const int idx = mfi.index();

            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                arr(i,j,k) = owner;
            });

            for (const auto& iv : pshifts)
            {
                ba.intersections(bx+iv, isects);                    
                for (const auto& is : isects)
                {
                    const int oi = is.first;
                    const Box& obx = is.second-iv;
                    if ((oi < idx) || (oi == idx && iv < IntVect::TheZeroVector())) 
                    {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(obx, i, j, k,
                        {
                            arr(i,j,k) = nonowner;
                        });
                    }
                }
            }
        }
    }

    return p;
}

}
