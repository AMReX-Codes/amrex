
#include <AMReX_BLassert.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_FabArrayUtility.H>
#include <AMReX_REAL.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#include <algorithm>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <map>
#include <limits>

namespace amrex {

namespace
{
    bool initialized = false;
#ifdef AMREX_MEM_PROFILING
    int num_multifabs     = 0;
    int num_multifabs_hwm = 0;
#endif
}

Real
MultiFab::Dot (const MultiFab& x, int xcomp,
               const MultiFab& y, int ycomp,
               int numcomp, int nghost, bool local)
{
    return amrex::Dot(x,xcomp,y,ycomp,numcomp,IntVect(nghost),local);
}

Real
MultiFab::Dot (const MultiFab& x, int xcomp, int numcomp, int nghost, bool local)
{
    BL_ASSERT(x.nGrow() >= nghost);

    BL_PROFILE("MultiFab::Dot()");

    Real sm = Real(0.0);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& xma = x.const_arrays();
        sm = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, x, IntVect(nghost),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            Real t = Real(0.0);
            auto const& xfab = xma[box_no];
            for (int n = 0; n < numcomp; ++n) {
                t += xfab(i,j,k,xcomp+n) * xfab(i,j,k,xcomp+n);
            }
            return t;
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(x,true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.growntilebox(nghost);
            Array4<Real const> const& xfab = x.const_array(mfi);
            AMREX_LOOP_4D(bx, numcomp, i, j, k, n,
            {
                sm += xfab(i,j,k,xcomp+n) * xfab(i,j,k,xcomp+n);
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
    }

    return sm;
}


Real
MultiFab::Dot (const iMultiFab& mask,
               const MultiFab& x, int xcomp,
               const MultiFab& y, int ycomp,
               int numcomp, int nghost, bool local)
{
    BL_ASSERT(x.boxArray() == y.boxArray());
    BL_ASSERT(x.boxArray() == mask.boxArray());
    BL_ASSERT(x.DistributionMap() == y.DistributionMap());
    BL_ASSERT(x.DistributionMap() == mask.DistributionMap());
    BL_ASSERT(x.nGrow() >= nghost && y.nGrow() >= nghost);
    BL_ASSERT(mask.nGrow() >= nghost);

    Real sm = Real(0.0);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& xma = x.const_arrays();
        auto const& yma = y.const_arrays();
        auto const& mma = mask.const_arrays();
        sm = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, x, IntVect(nghost),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            Real t = Real(0.0);
            if (mma[box_no](i,j,k)) {
                auto const& xfab = xma[box_no];
                auto const& yfab = yma[box_no];
                for (int n = 0; n < numcomp; ++n) {
                    t += xfab(i,j,k,xcomp+n) * yfab(i,j,k,ycomp+n);
                }
            }
            return t;
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(x,true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.growntilebox(nghost);
            Array4<Real const> const& xfab = x.const_array(mfi);
            Array4<Real const> const& yfab = y.const_array(mfi);
            Array4<int const> const& mfab = mask.const_array(mfi);
            AMREX_LOOP_4D(bx, numcomp, i, j, k, n,
            {
                if (mfab(i,j,k)) {
                    sm += xfab(i,j,k,xcomp+n) * yfab(i,j,k,ycomp+n);
                }
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
    }

    return sm;
}

void
MultiFab::Add (MultiFab& dst, const MultiFab& src,
               int srccomp, int dstcomp, int numcomp, int nghost)
{
    amrex::Add(dst,src,srccomp,dstcomp,numcomp,nghost);
}

void
MultiFab::Add (MultiFab& dst, const MultiFab& src,
               int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost) && src.nGrowVect().allGE(nghost));

    BL_PROFILE("MultiFab::Add()");

    amrex::Add(dst, src, srccomp, dstcomp, numcomp, nghost);
}

void
MultiFab::Copy (MultiFab& dst, const MultiFab& src,
                int srccomp, int dstcomp, int numcomp, int nghost)
{
    amrex::Copy(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::Copy (MultiFab& dst, const MultiFab& src,
                int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
// don't have to BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost));

    BL_PROFILE("MultiFab::Copy()");

    amrex::Copy(dst,src,srccomp,dstcomp,numcomp,nghost);
}

void
MultiFab::Swap (MultiFab& dst, MultiFab& src,
                int srccomp, int dstcomp, int numcomp, int nghost)
{
    Swap(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::Swap (MultiFab& dst, MultiFab& src,
                int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost) && src.nGrowVect().allGE(nghost));

    BL_PROFILE("MultiFab::Swap()");

    // We can take a shortcut and do a std::swap if we're swapping all of the data
    // and they are allocated in the same Arena.

    bool explicit_swap = true;

    if (srccomp == dstcomp && dstcomp == 0 && src.nComp() == dst.nComp() &&
        src.nGrowVect() == nghost && src.nGrowVect() == dst.nGrowVect() &&
        src.arena() == dst.arena() && src.hasEBFabFactory() == dst.hasEBFabFactory()) {
        explicit_swap = false;
    }

    if (!explicit_swap) {

        std::swap(dst, src);

    } else {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion() && dst.isFusingCandidate()) {
            auto const& dstma = dst.arrays();
            auto const& srcma = src.arrays();
            ParallelFor(dst, nghost, numcomp,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
            {
                const amrex::Real tmp = dstma[box_no](i,j,k,n+dstcomp);
                dstma[box_no](i,j,k,n+dstcomp) = srcma[box_no](i,j,k,n+srccomp);
                srcma[box_no](i,j,k,n+srccomp) = tmp;
            });
            if (!Gpu::inNoSyncRegion()) {
                Gpu::streamSynchronize();
            }
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox(nghost);
                if (bx.ok()) {
                    auto sfab = src.array(mfi);
                    auto dfab = dst.array(mfi);
                    AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, numcomp, i, j, k, n,
                    {
                        const amrex::Real tmp = dfab(i,j,k,n+dstcomp);
                        dfab(i,j,k,n+dstcomp) = sfab(i,j,k,n+srccomp);
                        sfab(i,j,k,n+srccomp) = tmp;
                    });
                }
            }
        }
    }
}


void
MultiFab::Subtract (MultiFab& dst, const MultiFab& src,
                    int srccomp, int dstcomp, int numcomp, int nghost)
{
    Subtract(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::Subtract (MultiFab& dst, const MultiFab& src,
                    int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost) && src.nGrowVect().allGE(nghost));

    BL_PROFILE("MultiFab::Subtract()");

    amrex::Subtract(dst,src,srccomp,dstcomp,numcomp,nghost);
}

void
MultiFab::Multiply (MultiFab& dst, const MultiFab& src,
                    int srccomp, int dstcomp, int numcomp, int nghost)
{
    Multiply(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::Multiply (MultiFab& dst, const MultiFab& src,
                    int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost) && src.nGrowVect().allGE(nghost));

    BL_PROFILE("MultiFab::Multiply()");

    amrex::Multiply(dst,src,srccomp,dstcomp,numcomp,nghost);
}

void
MultiFab::Divide (MultiFab& dst, const MultiFab& src,
                  int srccomp, int dstcomp, int numcomp, int nghost)
{
    Divide(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::Divide (MultiFab& dst, const MultiFab& src,
                  int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost) && src.nGrowVect().allGE(nghost));

    BL_PROFILE("MultiFab::Divide()");

    amrex::Divide(dst,src,srccomp,dstcomp,numcomp,nghost);
}

void
MultiFab::Saxpy (MultiFab& dst, Real a, const MultiFab& src,
                 int srccomp, int dstcomp, int numcomp, int nghost)
{
    Saxpy(dst,a,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::Xpay (MultiFab& dst, Real a, const MultiFab& src,
                int srccomp, int dstcomp, int numcomp, int nghost)
{
    Xpay(dst,a,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::LinComb (MultiFab& dst,
                   Real a, const MultiFab& x, int xcomp,
                   Real b, const MultiFab& y, int ycomp,
                   int dstcomp, int numcomp, int nghost)
{
    LinComb(dst,a,x,xcomp,b,y,ycomp,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::AddProduct (MultiFab& dst,
                      const MultiFab& src1, int comp1,
                      const MultiFab& src2, int comp2,
                      int dstcomp, int numcomp, int nghost)
{
    AddProduct(dst,src1,comp1,src2,comp2,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::AddProduct (MultiFab& dst,
                      const MultiFab& src1, int comp1,
                      const MultiFab& src2, int comp2,
                      int dstcomp, int numcomp, const IntVect& nghost)
{
    BL_ASSERT(dst.boxArray() == src1.boxArray());
    BL_ASSERT(dst.distributionMap == src1.distributionMap);
    BL_ASSERT(dst.boxArray() == src2.boxArray());
    BL_ASSERT(dst.distributionMap == src2.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost) && src1.nGrowVect().allGE(nghost) && src2.nGrowVect().allGE(nghost));

    BL_PROFILE("MultiFab::AddProduct()");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && dst.isFusingCandidate()) {
        auto const& dstma = dst.arrays();
        auto const& src1ma = src1.const_arrays();
        auto const& src2ma = src2.const_arrays();
        ParallelFor(dst, nghost, numcomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            dstma[box_no](i,j,k,n+dstcomp) += src1ma[box_no](i,j,k,n+comp1)
                *                             src2ma[box_no](i,j,k,n+comp2);
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            if (bx.ok()) {
                auto const s1fab = src1.array(mfi);
                auto const s2fab = src2.array(mfi);
                auto        dfab =  dst.array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, numcomp, i, j, k, n,
                {
                    dfab(i,j,k,n+dstcomp) += s1fab(i,j,k,n+comp1) * s2fab(i,j,k,n+comp2);
                });
            }
        }
    }
}

void
MultiFab::plus (Real val, int  nghost)
{
    plus(val,0,n_comp,nghost);
}

void
MultiFab::plus (Real val, const Box& region, int nghost)
{
    plus(val,region,0,n_comp,nghost);
}

void
MultiFab::mult (Real val, int nghost)
{
    mult(val,0,n_comp,nghost);
}

void
MultiFab::mult (Real val, const Box& region, int nghost)
{
    mult(val,region,0,n_comp,nghost);
}

void
MultiFab::invert (Real numerator, int nghost)
{
    invert(numerator,0,n_comp,nghost);
}

void
MultiFab::invert (Real numerator, const Box& region, int nghost)
{
    invert(numerator,region,0,n_comp,nghost);
}

void
MultiFab::negate (int nghost)
{
    negate(0,n_comp,nghost);
}

void
MultiFab::negate (const Box& region, int nghost)
{
    negate(region,0,n_comp,nghost);
}

void
MultiFab::Initialize ()
{
    if (initialized) return;
    initialized = true;

    amrex::ExecOnFinalize(MultiFab::Finalize);

#ifdef AMREX_MEM_PROFILING
    MemProfiler::add("MultiFab", std::function<MemProfiler::NBuildsInfo()>
                     ([] () -> MemProfiler::NBuildsInfo {
                         return {num_multifabs, num_multifabs_hwm};
                     }));
#endif
}

void
MultiFab::Finalize ()
{
    initialized = false;
}

MultiFab::MultiFab () noexcept // NOLINT(modernize-use-equals-default)
{
#ifdef AMREX_MEM_PROFILING
    ++num_multifabs;
    num_multifabs_hwm = std::max(num_multifabs_hwm, num_multifabs);
#endif
}

MultiFab::MultiFab (Arena* a) noexcept
    : FabArray<FArrayBox>(a)
{
#ifdef AMREX_MEM_PROFILING
    ++num_multifabs;
    num_multifabs_hwm = std::max(num_multifabs_hwm, num_multifabs);
#endif
}

MultiFab::MultiFab (const BoxArray&            bxs,
                    const DistributionMapping& dm,
                    int                        ncomp,
                    int                        ngrow,
                    const MFInfo&              info,
                    const FabFactory<FArrayBox>& factory)
    : MultiFab(bxs,dm,ncomp,IntVect(ngrow),info,factory)
{}

MultiFab::MultiFab (const BoxArray&            bxs,
                    const DistributionMapping& dm,
                    int                        ncomp,
                    const IntVect&             ngrow,
                    const MFInfo&              info,
                    const FabFactory<FArrayBox>& factory)
    :
    FabArray<FArrayBox>(bxs,dm,ncomp,ngrow,info,factory)
{
    if (SharedMemory() && info.alloc) initVal();  // else already done in FArrayBox
#ifdef AMREX_MEM_PROFILING
    ++num_multifabs;
    num_multifabs_hwm = std::max(num_multifabs_hwm, num_multifabs);
#endif
}

MultiFab::MultiFab (const MultiFab& rhs, MakeType maketype, int scomp, int ncomp)
    :
    FabArray<FArrayBox>(rhs, maketype, scomp, ncomp)
{
#ifdef AMREX_MEM_PROFILING
    ++num_multifabs;
    num_multifabs_hwm = std::max(num_multifabs_hwm, num_multifabs);
#endif
}

MultiFab::MultiFab (MultiFab&& rhs) noexcept
    : FabArray<FArrayBox>(std::move(rhs))
{
#ifdef AMREX_MEM_PROFILING
    ++num_multifabs;
    num_multifabs_hwm = std::max(num_multifabs_hwm, num_multifabs);
#endif
}

MultiFab::~MultiFab() // NOLINT(modernize-use-equals-default)
{
#ifdef AMREX_MEM_PROFILING
    --num_multifabs;
#endif
}

MultiFab&
MultiFab::operator= (Real r)
{
    setVal(r);
    return *this;
}

void
MultiFab::define (const BoxArray&            bxs,
                  const DistributionMapping& dm,
                  int                        nvar,
                  int                        ngrow,
                  const MFInfo&              info,
                  const FabFactory<FArrayBox>& factory)
{
    define(bxs, dm, nvar, IntVect(ngrow), info, factory);
    if (SharedMemory() && info.alloc) initVal();  // else already done in FArrayBox
}

void
MultiFab::define (const BoxArray&            bxs,
                  const DistributionMapping& dm,
                  int                        nvar,
                  const IntVect&             ngrow,
                  const MFInfo&              info,
                  const FabFactory<FArrayBox>& factory)
{
    this->FabArray<FArrayBox>::define(bxs,dm,nvar,ngrow,info,factory);
    if (SharedMemory() && info.alloc) initVal();  // else already done in FArrayBox
}

void
MultiFab::initVal ()
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*this)[mfi];
        fab.initVal();
    }
}

bool
MultiFab::contains_nan (int scomp, int ncomp, int ngrow, bool local) const
{
    return contains_nan(scomp, ncomp, IntVect(ngrow), local);
}

bool
MultiFab::contains_nan (int scomp, int ncomp, const IntVect& ngrow, bool local) const
{
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(scomp + ncomp <= nComp());
    BL_ASSERT(ncomp >  0 && ncomp <= nComp());
    BL_ASSERT(IntVect::TheZeroVector().allLE(ngrow) && ngrow.allLE(nGrowVect()));

    BL_PROFILE("MultiFab::contains_nan()");

    bool r = false;
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& ma = this->const_arrays();
        r = ParReduce(TypeList<ReduceOpLogicalOr>{}, TypeList<bool>{}, *this, ngrow, ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept -> GpuTuple<bool>
        {
            return amrex::isnan(ma[box_no](i,j,k,n+scomp));
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(||:r)
#endif
        for (MFIter mfi(*this,true); mfi.isValid() && !r; ++mfi)
        {
            Box const& bx = mfi.growntilebox(ngrow);
            Array4<Real const> const& fab = this->const_array(mfi);
            AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
            {
                r = r || amrex::isnan(fab(i,j,k,n+scomp));
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Or(r, ParallelContext::CommunicatorSub());
    }

    return r;
}

bool
MultiFab::contains_nan (bool local) const
{
    return contains_nan(0,nComp(),nGrowVect(),local);
}

bool
MultiFab::contains_inf (int scomp, int ncomp, IntVect const& ngrow, bool local) const
{
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(scomp + ncomp <= nComp());
    BL_ASSERT(ncomp >  0 && ncomp <= nComp());
    BL_ASSERT(IntVect::TheZeroVector().allLE(ngrow) && ngrow.allLE(nGrowVect()));

    BL_PROFILE("MultiFab::contains_inf()");

    bool r = false;
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& ma = this->const_arrays();
        r = ParReduce(TypeList<ReduceOpLogicalOr>{}, TypeList<bool>{}, *this, ngrow, ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept -> GpuTuple<bool>
        {
            return amrex::isinf(ma[box_no](i,j,k,n+scomp));
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(||:r)
#endif
        for (MFIter mfi(*this,true); mfi.isValid() && !r; ++mfi)
        {
            Box const& bx = mfi.growntilebox(ngrow);
            Array4<Real const> const& fab = this->const_array(mfi);
            AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
            {
                r = r || amrex::isinf(fab(i,j,k,n+scomp));
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Or(r, ParallelContext::CommunicatorSub());
    }

    return r;
}

bool
MultiFab::contains_inf (int scomp, int ncomp, int ngrow, bool local) const
{
    return contains_inf(scomp,ncomp,IntVect(ngrow),local);
}

bool
MultiFab::contains_inf (bool local) const
{
    return contains_inf(0,nComp(),nGrow(),local);
}

Real
MultiFab::min (int comp, int nghost, bool local) const
{
    BL_PROFILE("MultiFab::min()");

    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    Real mn = std::numeric_limits<Real>::max();

#ifdef AMREX_USE_EB
    if ( this->hasEBFabFactory() )
    {
        const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(this->Factory());
        auto const& flags = ebfactory.getMultiEBCellFlagFab();
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            auto const& flagsma = flags.const_arrays();
            auto const& ma = this->const_arrays();
            mn = ParReduce(TypeList<ReduceOpMin>{}, TypeList<Real>{}, *this, IntVect(nghost),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
            {
                if (flagsma[box_no](i,j,k).isCovered()) {
                    return AMREX_REAL_MAX;
                } else {
                    return ma[box_no](i,j,k,comp);
                }
            });
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(min:mn)
#endif
            for (MFIter mfi(*this,true); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.growntilebox(nghost);
                if (flags[mfi].getType(bx) != FabType::covered) {
                    auto const& flag = flags.const_array(mfi);
                    auto const& a = this->const_array(mfi);
                    AMREX_LOOP_3D(bx, i, j, k,
                    {
                        if (!flag(i,j,k).isCovered()) {
                            mn = std::min(mn, a(i,j,k,comp));
                        }
                    });
                }
            }
        }
    }
    else
#endif
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            auto const& ma = this->const_arrays();
            mn = ParReduce(TypeList<ReduceOpMin>{}, TypeList<Real>{}, *this, IntVect(nghost),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
            {
                return ma[box_no](i,j,k,comp);
            });
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(min:mn)
#endif
            for (MFIter mfi(*this,true); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.growntilebox(nghost);
                auto const& a = this->const_array(mfi);
                AMREX_LOOP_3D(bx, i, j, k,
                {
                    mn = std::min(mn, a(i,j,k,comp));
                });
            }
        }
    }

    if (!local) {
        ParallelAllReduce::Min(mn, ParallelContext::CommunicatorSub());
    }

    return mn;
}

Real
MultiFab::min (const Box& region, int comp, int nghost, bool local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    BL_PROFILE("MultiFab::min(region)");

    Real mn = std::numeric_limits<Real>::max();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& ma = this->const_arrays();
        mn = ParReduce(TypeList<ReduceOpMin>{}, TypeList<Real>{}, *this, IntVect(nghost),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            if (region.contains(i,j,k)) {
                return ma[box_no](i,j,k,comp);
            } else {
                return AMREX_REAL_MAX;
            }
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(min:mn)
#endif
        for (MFIter mfi(*this,true); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.growntilebox(nghost) & region;
            if (bx.ok()) {
                auto const& a = this->const_array(mfi);
                AMREX_LOOP_3D(bx, i, j, k,
                {
                    mn = std::min(mn, a(i,j,k,comp));
                });
            }
        }
    }

    if (!local) {
        ParallelAllReduce::Min(mn, ParallelContext::CommunicatorSub());
    }

    return mn;
}

Real
MultiFab::max (int comp, int nghost, bool local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    BL_PROFILE("MultiFab::max()");

    Real mx = std::numeric_limits<Real>::lowest();

#ifdef AMREX_USE_EB
    if ( this->hasEBFabFactory() )
    {
        const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(this->Factory());
        auto const& flags = ebfactory.getMultiEBCellFlagFab();
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            auto const& flagsma = flags.const_arrays();
            auto const& ma = this->const_arrays();
            mx = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, *this, IntVect(nghost),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
            {
                if (flagsma[box_no](i,j,k).isCovered()) {
                    return AMREX_REAL_LOWEST;
                } else {
                    return ma[box_no](i,j,k,comp);
                }
            });
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:mx)
#endif
            for (MFIter mfi(*this,true); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.growntilebox(nghost);
                if (flags[mfi].getType(bx) != FabType::covered) {
                    auto const& flag = flags.const_array(mfi);
                    auto const& a = this->const_array(mfi);
                    AMREX_LOOP_3D(bx, i, j, k,
                    {
                        if (!flag(i,j,k).isCovered()) {
                            mx = std::max(mx, a(i,j,k,comp));
                        }
                    });
                }
            }
        }
    }
    else
#endif
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            auto const& ma = this->const_arrays();
            mx = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, *this, IntVect(nghost),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
            {
                return ma[box_no](i,j,k,comp);
            });
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:mx)
#endif
            for (MFIter mfi(*this,true); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.growntilebox(nghost);
                auto const& a = this->const_array(mfi);
                AMREX_LOOP_3D(bx, i, j, k,
                {
                    mx = std::max(mx, a(i,j,k,comp));
                });
            }
        }
    }

    if (!local) {
        ParallelAllReduce::Max(mx, ParallelContext::CommunicatorSub());
    }

    return mx;
}

Real
MultiFab::max (const Box& region, int comp, int nghost, bool local) const
{
    BL_PROFILE("MultiFab::max(region)");

    Real mx = std::numeric_limits<Real>::lowest();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& ma = this->const_arrays();
        mx = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, *this, IntVect(nghost),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            if (region.contains(i,j,k)) {
                return ma[box_no](i,j,k,comp);
            } else {
                return AMREX_REAL_LOWEST;
            }
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:mx)
#endif
        for (MFIter mfi(*this,true); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.growntilebox(nghost) & region;
            if (bx.ok()) {
                auto const& a = this->const_array(mfi);
                AMREX_LOOP_3D(bx, i, j, k,
                {
                    mx = std::max(mx, a(i,j,k,comp));
                });
            }
        }
    }

    if (!local) {
        ParallelAllReduce::Max(mx, ParallelContext::CommunicatorSub());
    }

    return mx;
}

namespace {

IntVect
indexFromValue (MultiFab const& mf, int comp, int nghost, Real value, MPI_Op mmloc)
{
    IntVect loc = indexFromValue(mf, comp, IntVect{nghost}, value);

#ifdef BL_USE_MPI
    const int NProcs = ParallelContext::NProcsSub();
    if (NProcs > 1)
    {
        struct {
            Real mm;
            int rank;
        } in, out;
        in.mm = value;
        in.rank = ParallelContext::MyProcSub();
        MPI_Datatype datatype = (sizeof(Real) == sizeof(double))
            ? MPI_DOUBLE_INT : MPI_FLOAT_INT;
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
MultiFab::minIndex (int comp, int nghost) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    Real mn = this->min(comp, nghost, true);
    return indexFromValue(*this, comp, nghost, mn, MPI_MINLOC);
}

IntVect
MultiFab::maxIndex (int comp, int nghost) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    Real mx = this->max(comp, nghost, true);
    return indexFromValue(*this, comp, nghost, mx, MPI_MAXLOC);
}

Real
MultiFab::norm0 (const iMultiFab& mask, int comp, int nghost, bool local) const
{
    return FabArray<FArrayBox>::norminf(mask, comp, 1, IntVect(nghost), local);
}

Real
MultiFab::norm0 (int comp, int nghost, bool local, bool ignore_covered ) const
{
    return FabArray<FArrayBox>::norminf(comp, 1, IntVect(nghost), local, ignore_covered);
}

Real
MultiFab::norm0 (int comp, int ncomp, IntVect const& nghost, bool local, bool ignore_covered ) const
{
    return FabArray<FArrayBox>::norminf(comp, ncomp, nghost, local, ignore_covered);
}

Vector<Real>
MultiFab::norm0 (const Vector<int>& comps, int nghost, bool local, bool ignore_covered) const
{
    int n = static_cast<int>(comps.size());
    Vector<Real> nm0;
    nm0.reserve(n);

    for (int comp : comps) {
        nm0.push_back(this->norm0(comp, nghost, true, ignore_covered));
    }

    if (!local) {
        ParallelAllReduce::Max(nm0.dataPtr(), n, ParallelContext::CommunicatorSub());
    }

    return nm0;
}

Real
MultiFab::norm2 (int comp) const
{
    BL_ASSERT(ixType().cellCentered());

    Real nm2 = MultiFab::Dot(*this, comp, 1, 0);
    nm2 = std::sqrt(nm2);
    return nm2;
}

Real
MultiFab::norm2 (int comp, const Periodicity& period) const
{
    BL_PROFILE("MultiFab::norm2(period)");

    Real nm2 = Real(0.0);

    auto mask = OverlapMask(period);

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& ma = this->const_arrays();
        auto const& maskma = mask->const_arrays();
        nm2 = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, *this,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            Real tmp = ma[box_no](i,j,k,comp);
            return tmp*tmp/maskma[box_no](i,j,k);
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:nm2)
#endif
        for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            auto const& a = this->const_array(mfi);
            auto const& m = mask->const_array(mfi);
            AMREX_LOOP_3D(bx, i, j, k,
            {
                Real tmp = a(i,j,k,comp);
                nm2 += tmp*tmp/m(i,j,k);
            });
        }
    }

    ParallelAllReduce::Sum(nm2, ParallelContext::CommunicatorSub());
    return std::sqrt(nm2);
}

Vector<Real>
MultiFab::norm2 (const Vector<int>& comps) const
{
    BL_ASSERT(ixType().cellCentered());

    int n = static_cast<int>(comps.size());
    Vector<Real> nm2;
    nm2.reserve(n);

    for (int comp : comps) {
        nm2.push_back(this->norm2(comp));
    }

    return nm2;
}

Real
MultiFab::norm1 (int comp, const Periodicity& period, bool ignore_covered ) const
{
    amrex::ignore_unused(ignore_covered);

    MultiFab tmpmf(this->boxArray(), this->DistributionMap(), 1, 0,
                   MFInfo(), this->Factory());

    MultiFab::Copy(tmpmf, *this, comp, 0, 1, 0);

#ifdef AMREX_USE_EB
    if ( this -> hasEBFabFactory() && ignore_covered )
        EB_set_covered( tmpmf, Real(0.0) );
#endif

    auto mask = OverlapMask(period);
    MultiFab::Divide(tmpmf, *mask, 0, 0, 1, 0);

    return tmpmf.norm1(0, 0);
}

Real
MultiFab::norm1 (int comp, int ngrow, bool local) const
{
    BL_PROFILE("MultiFab::norm1");

    Real nm1 = Real(0.0);

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& ma = this->const_arrays();
        nm1 = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, *this, IntVect(ngrow),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            return std::abs(ma[box_no](i,j,k,comp));
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(+:nm1)
#endif
        for (MFIter mfi(*this,true); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.growntilebox(ngrow);
            auto const& a = this->const_array(mfi);
            AMREX_LOOP_3D(bx, i, j, k,
            {
                nm1 += std::abs(a(i,j,k,comp));
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Sum(nm1, ParallelContext::CommunicatorSub());
    }

    return nm1;
}

Vector<Real>
MultiFab::norm1 (const Vector<int>& comps, int ngrow, bool local) const
{
    BL_ASSERT(ixType().cellCentered());

    int n = static_cast<int>(comps.size());
    Vector<Real> nm1;
    nm1.reserve(n);

    for (int comp : comps) {
        nm1.push_back(this->norm1(comp, ngrow, true));
    }

    if (!local)
        ParallelAllReduce::Sum(nm1.dataPtr(), n, ParallelContext::CommunicatorSub());

    return nm1;
}

Real
MultiFab::sum (int comp, bool local) const
{
    return FabArray<FArrayBox>::sum(comp, IntVect(0), local);
}

Real
MultiFab::sum_unique (int comp,
                      bool local,
                      const Periodicity& period) const
{
    BL_PROFILE("MultiFab::sum_unique()");

    // no duplicatly distributed points if cell centered
    if (ixType().cellCentered())
        return this->sum(comp, local);

    // Owner is the grid with the lowest grid number containing the data
    std::unique_ptr<iMultiFab> owner_mask = OwnerMask(period);

    Real sm = Real(0.0);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& ma = this->const_arrays();
        auto const& msk = owner_mask->const_arrays();
        sm = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, *this, IntVect(0),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                       -> GpuTuple<Real>
        {
            return msk[box_no](i,j,k) ? ma[box_no](i,j,k,comp) : 0.0_rt;
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real const> const& a = this->const_array(mfi);
            Array4<int const> const& msk = owner_mask->const_array(mfi);
            Real tmp = 0.0_rt;
            AMREX_LOOP_3D(bx, i, j, k,
            {
                tmp += msk(i,j,k) ? a(i,j,k,comp) : 0.0_rt;
            });
            sm += tmp; // Do it this way so that it does not break regression tests.
        }
    }

    if (!local) {
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
    }

    return sm;
}

void
MultiFab::minus (const MultiFab& mf, int strt_comp, int num_comp, int nghost)
{
    MultiFab::Subtract(*this, mf, strt_comp, strt_comp, num_comp, nghost);
}

void
MultiFab::divide (const MultiFab& mf, int strt_comp, int num_comp, int nghost)
{
    MultiFab::Divide(*this, mf, strt_comp, strt_comp, num_comp, nghost);
}

void
MultiFab::plus (Real val, int comp, int num_comp, int nghost)
{

    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    FabArray<FArrayBox>::plus(val, comp, num_comp, nghost);
}

void
MultiFab::plus (Real val, const Box& region, int comp, int num_comp, int nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    FabArray<FArrayBox>::plus(val,region,comp,num_comp,nghost);
}

void
MultiFab::plus (const MultiFab& mf, int strt_comp, int num_comp, int nghost)
{
    MultiFab::Add(*this, mf, strt_comp, strt_comp, num_comp, nghost);
}

void
MultiFab::mult (Real val, int comp, int num_comp, int  nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    FabArray<FArrayBox>::mult(val,comp,num_comp,nghost);
}

void
MultiFab::mult (Real val, const Box& region, int comp, int num_comp, int nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    FabArray<FArrayBox>::mult(val,region,comp,num_comp,nghost);
}

void
MultiFab::invert (Real numerator, int comp, int num_comp, int nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    FabArray<FArrayBox>::invert(numerator,comp,num_comp,nghost);
}

void
MultiFab::invert (Real numerator, const Box& region, int comp, int num_comp, int nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

    FabArray<FArrayBox>::invert(numerator,region,comp,num_comp,nghost);
}

void
MultiFab::negate (int comp, int num_comp, int nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);

    FabArray<FArrayBox>::mult(-1., comp, num_comp, nghost);
}

void
MultiFab::negate (const Box& region, int comp, int num_comp, int nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);

    FabArray<FArrayBox>::mult(-1.,region,comp,num_comp,nghost);
}

std::unique_ptr<MultiFab>
MultiFab::OverlapMask (const Periodicity& period) const
{
    BL_PROFILE("MultiFab::OverlapMask()");

    const BoxArray& ba = boxArray();
    const DistributionMapping& dm = DistributionMap();

    auto p = std::make_unique<MultiFab>(ba,dm,1,0, MFInfo(), Factory());

    const std::vector<IntVect>& pshifts = period.shiftIntVect();

    Vector<Array4BoxTag<Real> > tags;

    bool run_on_gpu = Gpu::inLaunchRegion();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!run_on_gpu)
#endif
    {
        std::vector< std::pair<int,Box> > isects;

        for (MFIter mfi(*p); mfi.isValid(); ++mfi)
        {
            const Box& bx = (*p)[mfi].box();
            Array4<Real> const& arr = p->array(mfi);

            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                arr(i,j,k) = Real(0.0);
            });

            for (const auto& iv : pshifts)
            {
                ba.intersections(bx+iv, isects);
                for (const auto& is : isects)
                {
                    Box const& b = is.second-iv;
                    if (run_on_gpu) {
                        tags.push_back({arr,b});
                    } else {
                        amrex::LoopConcurrentOnCpu(b, [=] (int i, int j, int k) noexcept
                        {
                            arr(i,j,k) += Real(1.0);
                        });
                    }
                }
            }
        }
    }

#ifdef AMREX_USE_GPU
    amrex::ParallelFor(tags, 1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n, Array4BoxTag<Real> const& tag) noexcept
    {
        Real* p = tag.dfab.ptr(i,j,k,n);
        Gpu::Atomic::AddNoRet(p, Real(1.0));
    });
#endif

    return p;
}

std::unique_ptr<iMultiFab>
MultiFab::OwnerMask (const Periodicity& period) const
{
    return amrex::OwnerMask(*this, period);
}

void
MultiFab::AverageSync (const Periodicity& period)
{
    BL_PROFILE("MultiFab::AverageSync()");

    if (ixType().cellCentered()) return;
    auto wgt = this->OverlapMask(period);
    wgt->invert(1.0, 0, 1);
    this->WeightedSync(*wgt, period);
}

void
MultiFab::WeightedSync (const MultiFab& wgt, const Periodicity& period)
{
    BL_PROFILE("MultiFab::WeightedSync()");

    if (ixType().cellCentered()) return;

    const int ncomp = nComp();
    for (int comp = 0; comp < ncomp; ++comp)
    {
        MultiFab::Multiply(*this, wgt, 0, comp, 1, 0);
    }

    MultiFab tmpmf(boxArray(), DistributionMap(), ncomp, 0, MFInfo(), Factory());
    tmpmf.setVal(Real(0.0));
    tmpmf.ParallelCopy(*this, period, FabArrayBase::ADD);

    MultiFab::Copy(*this, tmpmf, 0, 0, ncomp, 0);
}

void
MultiFab::OverrideSync (const iMultiFab& msk, const Periodicity& period)
{
    amrex::OverrideSync(*this, msk, period);
}

}
