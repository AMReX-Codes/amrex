
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
iMultiFab::plus (int val, int nghost)
{
    plus(val,0,n_comp,nghost);
}

void
iMultiFab::plus (int val, const Box& region, int nghost)
{
    plus(val,region,0,n_comp,nghost);
}

void
iMultiFab::mult (int val, int nghost)
{
    mult(val,0,n_comp,nghost);
}

void
iMultiFab::mult (int val, const Box& region, int nghost)
{
    mult(val,region,0,n_comp,nghost);
}

void
iMultiFab::negate (int nghost)
{
    negate(0,n_comp,nghost);
}

void
iMultiFab::negate (const Box& region, int nghost)
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

iMultiFab::~iMultiFab()
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
iMultiFab::min (int comp, int nghost, bool local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    int mn = amrex::ReduceMin(*this, nghost,
    [=] AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<int const> const& fab) -> int
    {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
        int r = std::numeric_limits<int>::max();
#else
        int r = INT_MAX;
#endif
        AMREX_LOOP_3D(bx, i, j, k,
        {
            r = amrex::min(r, fab(i,j,k,comp));
        });
        return r;
    });                          

    if (!local)
	ParallelDescriptor::ReduceIntMin(mn);

    return mn;
}

int
iMultiFab::min (const Box& region, int comp, int nghost, bool local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    int mn = amrex::ReduceMin(*this, nghost,
    [=] AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<int const> const& fab) -> int
    {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
        int r = std::numeric_limits<int>::max();
#else
        int r = INT_MAX;
#endif
        const Box& b = bx & region;
        AMREX_LOOP_3D(b, i, j, k,
        {
            r = amrex::min(r, fab(i,j,k,comp));
        });
        return r;
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
    [=] AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<int const> const& fab) -> int
    {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
        int r = std::numeric_limits<int>::lowest();
#else
        int r = INT_MIN;
#endif
        AMREX_LOOP_3D(bx, i, j, k,
        {
            r = amrex::max(r, fab(i,j,k,comp));
        });
        return r;
    });                          

    if (!local)
	ParallelDescriptor::ReduceIntMax(mx);

    return mx;
}

int
iMultiFab::max (const Box& region, int comp, int nghost, bool local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    int mx = amrex::ReduceMax(*this, nghost,
    [=] AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<int const> const& fab) -> int
    {
        const Box& b = bx & region;
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
        int r = std::numeric_limits<int>::lowest();
#else
        int r = INT_MIN;
#endif
        AMREX_LOOP_3D(b, i, j, k,
        {
            r = amrex::max(r, fab(i,j,k,comp));
        });
        return r;
    });

    if (!local)
	ParallelDescriptor::ReduceIntMax(mx);

    return mx;
}

Long
iMultiFab::sum (int comp, int nghost, bool local) const
{
    AMREX_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    Long sm = 0;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Long> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(*this); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            const auto& arr = this->const_array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return { static_cast<Long>(arr(i,j,k,comp)) };
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
            Array4<int const> const& fab = this->const_array(mfi);
            AMREX_LOOP_3D(bx, i, j, k,
            {
                sm += fab(i,j,k,comp);
            });
        }
    }

    if (!local) ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());

    return sm;
}

namespace {

static IntVect
indexFromValue (iMultiFab const& mf, int comp, int nghost, int value, MPI_Op mmloc)
{
    IntVect loc = indexFromValue(mf, comp, IntVect{nghost}, value);

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
#else
    amrex::ignore_unused(mmloc);
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
OwnerMask (FabArrayBase const& mf, const Periodicity& period, const IntVect& ngrow)
{
    BL_PROFILE("OwnerMask()");

    const BoxArray& ba = mf.boxArray();
    const DistributionMapping& dm = mf.DistributionMap();

    const int owner = 1;
    const int nonowner = 0;

    std::unique_ptr<iMultiFab> p{new iMultiFab(ba,dm,1,ngrow, MFInfo(),
                                               DefaultFabFactory<IArrayBox>())};
    const std::vector<IntVect>& pshifts = period.shiftIntVect();

    Vector<Array4BoxTag<int> > tags;

    bool run_on_gpu = Gpu::inLaunchRegion();
#ifdef _OPENMP
#pragma omp parallel if (!run_on_gpu)
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
                ba.intersections(bx+iv, isects, false, ngrow);
                for (const auto& is : isects)
                {
                    const int oi = is.first;
                    const Box& obx = is.second-iv;
                    if ((oi < idx) || (oi == idx && iv < IntVect::TheZeroVector())) 
                    {
                        if (run_on_gpu) {
                            tags.push_back({arr,obx});
                        } else {
                            // cannot use amrex::Loop because of a gcc bug.
                            const auto lo = amrex::lbound(obx);
                            const auto hi = amrex::ubound(obx);
                            for (int k = lo.z; k <= hi.z; ++k) {
                            for (int j = lo.y; j <= hi.y; ++j) {
                            AMREX_PRAGMA_SIMD
                            for (int i = lo.x; i <= hi.x; ++i) {
                                arr(i,j,k) = nonowner;
                            }}}
                        }
                    }
                }
            }
        }
    }

#ifdef AMREX_USE_GPU
    amrex::ParallelFor(tags, 1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n, Array4<int> const& a) noexcept
    {
        a(i,j,k,n) = nonowner;
    });
#endif

    return p;
}

}
