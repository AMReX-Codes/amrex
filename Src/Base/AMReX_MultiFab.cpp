
#include <algorithm>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <map>
#include <limits>

#include <AMReX_BLassert.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_iMultiFab.H>

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

namespace amrex {

namespace
{
    bool initialized = false;
#ifdef BL_MEM_PROFILING
    int num_multifabs     = 0;
    int num_multifabs_hwm = 0;
#endif
}

Real
MultiFab::Dot (const MultiFab& x, int xcomp,
	       const MultiFab& y, int ycomp,
	       int numcomp, int nghost, bool local)
{
    BL_ASSERT(x.boxArray() == y.boxArray());
    BL_ASSERT(x.DistributionMap() == y.DistributionMap());
    BL_ASSERT(x.nGrow() >= nghost && y.nGrow() >= nghost);

    BL_PROFILE("MultiFab::Dot()");

    Real sm = 0.0;

#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction && Gpu::notInLaunchRegion()) reduction(+:sm)
#endif
    {
        amrex::Gpu::DeviceScalar<Real> local_sm(0.0);
        Real* p = local_sm.dataPtr();
        for (MFIter mfi(x,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            FArrayBox const* xfab = &(x[mfi]);
            FArrayBox const* yfab = &(y[mfi]);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
            {
                Real t = xfab->dot(tbx, xcomp, *yfab, tbx, ycomp, numcomp);
                amrex::Gpu::Atomic::Add(p, t);
            });
        }
        sm += local_sm.dataValue();
    }

    if (!local) ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());

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

    Real sm = 0.0;

#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction && Gpu::notInLaunchRegion()) reduction(+:sm)
#endif
    {
        amrex::Gpu::DeviceScalar<Real> local_sm(0.0);
        Real* p = local_sm.dataPtr();
        for (MFIter mfi(x,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            FArrayBox const* xfab = &(x[mfi]);
            FArrayBox const* yfab = &(y[mfi]);
            IArrayBox const* mskfab = &(mask[mfi]);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
            {
                Real t = xfab->dotmask(*mskfab, tbx, xcomp, *yfab, tbx, ycomp, numcomp);
                amrex::Gpu::Atomic::Add(p, t);
            });
        }
        sm += local_sm.dataValue();
    }

    if (!local)
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());

    return sm;
}

void
MultiFab::Add (MultiFab&       dst,
	       const MultiFab& src,
	       int             srccomp,
	       int             dstcomp,
	       int             numcomp,
	       int             nghost)
{
    Add(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::Add (MultiFab&       dst,
	       const MultiFab& src,
	       int             srccomp,
	       int             dstcomp,
	       int             numcomp,
	       const IntVect&  nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost) && src.nGrowVect().allGE(nghost));

    BL_PROFILE("MultiFab::Add()");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        FArrayBox const* srcFab = &(src[mfi]);
        FArrayBox      * dstFab = &(dst[mfi]);

        if (bx.ok())
        {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
            {
                dstFab->plus(*srcFab, tbx, tbx, srccomp, dstcomp, numcomp);
            });
        }
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
    Copy(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::Copy (MultiFab&       dst,
                const MultiFab& src,
                int             srccomp,
                int             dstcomp,
                int             numcomp,
                const IntVect&  nghost)
{
// don't have to BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost));

    BL_PROFILE("MultiFab::Copy()");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        FArrayBox const* srcFab = &(src[mfi]);
        FArrayBox      * dstFab = &(dst[mfi]);

        if (bx.ok())
        {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
            {
                dstFab->copy(*srcFab, tbx, srccomp, tbx, dstcomp, numcomp);
            });
        }
    }
}


#ifdef USE_PERILLA
void
MultiFab::Copy (MultiFab&       dst,
                const MultiFab& src,
                int             f,
                int             srccomp,
                int             dstcomp,
                int             numcomp,
                const Box&      bx)
{
// don't have to    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    //BL_ASSERT(dst.nGrow() >= nghost); // && src.nGrow() >= nghost);

    int fis = src.IndexArray()[f];
    int fid = dst.IndexArray()[f];
    //const Box& bx = BoxLib::grow(dst[f].box(),nghost);
    //const Box& bx = dst[fid].box();

    if (bx.ok())
      dst[fid].copy(src[fid], bx, srccomp, bx, dstcomp, numcomp);

}
#endif


void
MultiFab::Subtract (MultiFab&       dst,
		    const MultiFab& src,
		    int             srccomp,
		    int             dstcomp,
		    int             numcomp,
		    int             nghost)
{
    Subtract(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::Subtract (MultiFab&       dst,
		    const MultiFab& src,
		    int             srccomp,
		    int             dstcomp,
		    int             numcomp,
		    const IntVect&  nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost) && src.nGrowVect().allGE(nghost));

    BL_PROFILE("MultiFab::Subtract()");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        FArrayBox const* srcFab = &(src[mfi]);
        FArrayBox      * dstFab = &(dst[mfi]);

        if (bx.ok())
        {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
            {
                dstFab->minus(*srcFab, tbx, tbx, srccomp, dstcomp, numcomp);
            });
        }
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
    Multiply(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::Multiply (MultiFab&       dst,
		    const MultiFab& src,
		    int             srccomp,
		    int             dstcomp,
		    int             numcomp,
		    const IntVect&  nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost) && src.nGrowVect().allGE(nghost));

    BL_PROFILE("MultiFab::Multiply()");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        FArrayBox const* srcFab = &(src[mfi]);
        FArrayBox      * dstFab = &(dst[mfi]);

        if (bx.ok())
        {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
            {
                dstFab->mult(*srcFab, tbx, tbx, srccomp, dstcomp, numcomp);
            });
        }
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
    Divide(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

void
MultiFab::Divide (MultiFab&       dst,
		  const MultiFab& src,
		  int             srccomp,
		  int             dstcomp,
		  int             numcomp,
		  const IntVect&  nghost)
{
    BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost) && src.nGrowVect().allGE(nghost));

    BL_PROFILE("MultiFab::Divide()");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        FArrayBox const* srcFab = &(src[mfi]);
        FArrayBox      * dstFab = &(dst[mfi]);

        if (bx.ok())
        {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
            {
                dstFab->divide(*srcFab, tbx, tbx, srccomp, dstcomp, numcomp);
            });
        }
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

    BL_PROFILE("MultiFab::Saxpy()");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);

        if (bx.ok()) {
            FArrayBox const* sfab = &(src[mfi]);
            FArrayBox      * dfab = &(dst[mfi]);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
            {
                dfab->saxpy(a, *sfab, tbx, tbx, srccomp, dstcomp, numcomp);
            });
        }
    }
}

void
MultiFab::Xpay (MultiFab&       dst,
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

    BL_PROFILE("MultiFab::Xpay()");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);

        if (bx.ok()) {
            FArrayBox const* sfab = &(src[mfi]);
            FArrayBox      * dfab = &(dst[mfi]);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
            {
                dfab->xpay(a, *sfab, tbx, tbx, srccomp, dstcomp, numcomp);
            });
        }
    }
}

void
MultiFab::LinComb (MultiFab&       dst,
		   Real            a,
		   const MultiFab& x,
		   int             xcomp,
		   Real            b,
		   const MultiFab& y,
		   int             ycomp,
		   int             dstcomp,
		   int             numcomp,
		   int             nghost)
{
    BL_ASSERT(dst.boxArray() == x.boxArray());
    BL_ASSERT(dst.distributionMap == x.distributionMap);
    BL_ASSERT(dst.boxArray() == y.boxArray());
    BL_ASSERT(dst.distributionMap == y.distributionMap);
    BL_ASSERT(dst.nGrow() >= nghost && x.nGrow() >= nghost && y.nGrow() >= nghost);

    BL_PROFILE("MultiFab::LinComb()");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
	
        if (bx.ok()) {
            FArrayBox const* xfab = &(x  [mfi]);
            FArrayBox const* yfab = &(y  [mfi]);
            FArrayBox      * dfab = &(dst[mfi]);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
            {
                dfab->linComb(*xfab,tbx,xcomp,*yfab,tbx,ycomp,a,b,tbx,dstcomp,numcomp);
            });
        }
    }
}

void
MultiFab::AddProduct (MultiFab&       dst,
		      const MultiFab& src1,
		      int             comp1,
		      const MultiFab& src2,
		      int             comp2,
		      int             dstcomp,
		      int             numcomp,
		      int             nghost)
{
    BL_ASSERT(dst.boxArray() == src1.boxArray());
    BL_ASSERT(dst.distributionMap == src1.distributionMap);
    BL_ASSERT(dst.boxArray() == src2.boxArray());
    BL_ASSERT(dst.distributionMap == src2.distributionMap);
    BL_ASSERT(dst.nGrow() >= nghost && src1.nGrow() >= nghost && src2.nGrow() >= nghost);

    BL_PROFILE("MultiFab::AddProduct()");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
	
        if (bx.ok()) {
            FArrayBox const* s1fab = &(src1[mfi]);
            FArrayBox const* s2fab = &(src2[mfi]);
            FArrayBox      *  dfab = &(dst [mfi]);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
            {
                dfab->addproduct(tbx, dstcomp, numcomp, *s1fab, comp1, *s2fab, comp2);
            });
        }
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
    initialized = true;

    amrex::ExecOnFinalize(MultiFab::Finalize);

#ifdef BL_MEM_PROFILING
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

MultiFab::MultiFab ()
{
#ifdef BL_MEM_PROFILING
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
#ifdef BL_MEM_PROFILING
    ++num_multifabs;
    num_multifabs_hwm = std::max(num_multifabs_hwm, num_multifabs);
#endif
}

MultiFab::MultiFab (const MultiFab& rhs, MakeType maketype, int scomp, int ncomp)
    :
    FabArray<FArrayBox>(rhs, maketype, scomp, ncomp)
{
#ifdef BL_MEM_PROFILING
    ++num_multifabs;
    num_multifabs_hwm = std::max(num_multifabs_hwm, num_multifabs);
#endif
}

MultiFab::MultiFab (const BoxArray& ba, const DistributionMapping& dm, int ncomp, int ngrow,
                    const Vector<Real*>& p)
    : MultiFab(ba, dm, ncomp, IntVect(ngrow), p)
{}

MultiFab::MultiFab (const BoxArray& ba, const DistributionMapping& dm, int ncomp, const IntVect& ngrow,
                    const Vector<Real*>& p)
    :
    FabArray<FArrayBox>(ba, dm, ncomp, ngrow, p)
{
#ifdef BL_MEM_PROFILING
    ++num_multifabs;
    num_multifabs_hwm = std::max(num_multifabs_hwm, num_multifabs);
#endif
}

MultiFab::MultiFab (MultiFab&& rhs) noexcept
    : FabArray<FArrayBox>(std::move(rhs))
{
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
    // Done in FabArray. Just Cuda wrapping and Tiling check here.
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this, TilingIfNotGPU()); mfi.isValid(); ++mfi)
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
    // TODO GPU -- CHECK
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(scomp + ncomp <= nComp());
    BL_ASSERT(ncomp >  0 && ncomp <= nComp());
    BL_ASSERT(ngrow >= 0 && ngrow <= nGrow());

    bool r = false;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(|:r)
#endif
    {
        amrex::Gpu::DeviceScalar<int> local_sm(0);
        int* p = local_sm.dataPtr();
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            if (local_sm.dataValue()) // Reduce calc if prevelant throughout MF.
            {
                const Box& bx = mfi.growntilebox(ngrow);
                FArrayBox const* fab = &(get(mfi));

                AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
                {	
                    int t = fab->contains_nan(tbx,scomp,ncomp);
                    amrex::Gpu::Atomic::Or(p, t);
                });
            }
        }
        r = (r || local_sm.dataValue());  // <---- Int to bool.
    }

    if (!local)
	ParallelAllReduce::Or(r, ParallelContext::CommunicatorSub());

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
    // TODO GPU -- CHECK
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(scomp + ncomp <= nComp());
    BL_ASSERT(ncomp >  0 && ncomp <= nComp());
    BL_ASSERT(ngrow >= 0 && ngrow <= nGrow());

    bool r = false;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(|:r)
#endif
    {
        amrex::Gpu::DeviceScalar<int> local_sm(0);
        int* p = local_sm.dataPtr();
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            if (local_sm.dataValue())  // Reduce calc if prevelant throughout MF.
            {
                const Box& bx = mfi.growntilebox(ngrow);
                FArrayBox const* fab = &(get(mfi));

                AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
                {	
                    int t = fab->contains_inf(tbx,scomp,ncomp);
                    amrex::Gpu::Atomic::Or(p, t);
                });
            }
        }
        r = (r || local_sm.dataValue());  // <---- Int to bool.
    }

    if (!local)
	ParallelAllReduce::Or(r, ParallelContext::CommunicatorSub());

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

bool 
MultiFab::is_nodal (int dir) const
{
    return boxArray().ixType().nodeCentered(dir);
}

bool 
MultiFab::is_cell_centered () const
{
    return boxArray().ixType().cellCentered();
}

Real
MultiFab::min (int comp,
               int nghost,
	       bool local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    Real mn = std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(min:mn) 
#endif
    {
        amrex::Gpu::DeviceScalar<Real> local_mn(std::numeric_limits<Real>::max());
        Real* p = local_mn.dataPtr();
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            FArrayBox const* fab = &(get(mfi));
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
            {
                Real t = fab->min(tbx, comp);
                amrex::Gpu::Atomic::Min(p, t);
            });
        }
        mn = std::min(mn, local_mn.dataValue());
    }

    if (!local)
	ParallelAllReduce::Min(mn, ParallelContext::CommunicatorSub());

    return mn;
}

Real
MultiFab::min (const Box& region,
               int        comp,
               int        nghost,
	       bool       local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    Real mn = std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(min:mn) 
#endif
    {
        amrex::Gpu::DeviceScalar<Real> local_mn(std::numeric_limits<Real>::max());
        Real* p = local_mn.dataPtr();
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost) & region;
            FArrayBox const* fab = &(get(mfi));
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
            {
                Real t = fab->min(tbx, comp);
                amrex::Gpu::Atomic::Min(p, t);
            });
        }
        mn = std::min(mn, local_mn.dataValue());
    }

    if (!local)
	ParallelAllReduce::Min(mn, ParallelContext::CommunicatorSub());

    return mn;

}

Real
MultiFab::max (int comp,
               int nghost,
	       bool local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    Real mx = -std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(max:mx) 
#endif
    {
        amrex::Gpu::DeviceScalar<Real> local_mx(-std::numeric_limits<Real>::max());
        Real* p = local_mx.dataPtr();
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            FArrayBox const* fab = &(get(mfi));
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
            {
                Real t = fab->max(tbx, comp);
                amrex::Gpu::Atomic::Max(p, t);
            });
        }
        mx = std::max(mx, local_mx.dataValue());
    }

    if (!local)
	ParallelAllReduce::Max(mx, ParallelContext::CommunicatorSub());

    return mx;
}

Real
MultiFab::max (const Box& region,
               int        comp,
               int        nghost,
	       bool       local) const
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    Real mx = -std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(max:mx) 
#endif
    {
        amrex::Gpu::DeviceScalar<Real> local_mx(-std::numeric_limits<Real>::max());
        Real* p = local_mx.dataPtr();
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost) & region;
            FArrayBox const* fab = &(get(mfi));
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
            {
                Real t = fab->max(tbx, comp);
                amrex::Gpu::Atomic::Max(p, t);
            });
        }
        mx = std::max(mx, local_mx.dataValue());
    }

    if (!local)
	ParallelAllReduce::Max(mx, ParallelContext::CommunicatorSub());

    return mx;
}

IntVect
MultiFab::minIndex (int comp,
                    int nghost) const
{
    // TODO GPU -- CHECK

    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    Real mn = std::numeric_limits<Real>::max();
    IntVect loc;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Real priv_mn = std::numeric_limits<Real>::max();
        IntVect priv_loc;

        amrex::Gpu::DeviceScalar<Real> local_mn(std::numeric_limits<Real>::max());
        Real* p = local_mn.dataPtr();
	for (MFIter mfi(*this); mfi.isValid(); ++mfi)
	{
	    const Box& bx = amrex::grow(mfi.validbox(),nghost);
            const FArrayBox* fab = &(get(mfi));

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
            {
                Real t = fab->min(tbx,comp);
                amrex::Gpu::Atomic::Min(p, t);
            });
	}
        priv_mn = std::min(priv_mn, local_mn.dataValue());
       

        amrex::Gpu::DeviceScalar<IntVect> local_loc(IntVect::TheZeroVector());
        IntVect* l = local_loc.dataPtr();
	for (MFIter mfi(*this); mfi.isValid(); ++mfi)
	{
	    const Box& bx = amrex::grow(mfi.validbox(),nghost);
            const FArrayBox* fab = &(get(mfi));

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
            {
                IntVect t_loc = fab->indexFromValue(priv_mn, tbx,comp);

                if (tbx.contains(t_loc))
                {
                    *l = t_loc; // For total safety, this should be a Gpu::Atomic.
                };
	    });
	}
        priv_loc = local_loc.dataValue();

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

    const int NProcs = ParallelContext::NProcsSub();
    if (NProcs > 1)
    {
        Vector<Real> mns(NProcs);
        Vector<int>  locs(NProcs * AMREX_SPACEDIM);

        auto comm = ParallelContext::CommunicatorSub();
        ParallelAllGather::AllGather(mn, mns.dataPtr(), comm);
        BL_ASSERT(sizeof(IntVect) == sizeof(int)*AMREX_SPACEDIM);
        ParallelAllGather::AllGather(loc.getVect(), AMREX_SPACEDIM, locs.dataPtr(), comm);

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

    return loc;
}

IntVect
MultiFab::maxIndex (int comp,
                    int nghost) const
{
    // TODO GPU -- CHECK

    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());

    Real mx = std::numeric_limits<Real>::lowest();
    IntVect loc;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Real priv_mx = std::numeric_limits<Real>::lowest();
        IntVect priv_loc;

        amrex::Gpu::DeviceScalar<Real> local_mx(std::numeric_limits<Real>::lowest());
        Real* p = local_mx.dataPtr();
	for (MFIter mfi(*this); mfi.isValid(); ++mfi)
	{
	    const Box& bx = amrex::grow(mfi.validbox(),nghost);
            const FArrayBox* fab = &(get(mfi));

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
            {
                Real t = fab->max(tbx,comp);
                amrex::Gpu::Atomic::Max(p, t);
            });
	}
        priv_mx = std::max(priv_mx, local_mx.dataValue());
       

        amrex::Gpu::DeviceScalar<IntVect> local_loc(IntVect::TheZeroVector());
        IntVect* l = local_loc.dataPtr();
	for (MFIter mfi(*this); mfi.isValid(); ++mfi)
	{
	    const Box& bx = amrex::grow(mfi.validbox(),nghost);
            const FArrayBox* fab = &(get(mfi));

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
            {
                IntVect t_loc = fab->indexFromValue(priv_mx, tbx,comp);

                if (tbx.contains(t_loc))
                {
                    *l = t_loc; // For total safety, this should be a Gpu::Atomic.
                };
	    });
	}
        priv_loc = local_loc.dataValue();

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

    const int NProcs = ParallelContext::NProcsSub();
    if (NProcs > 1)
    {
        Vector<Real> mxs(NProcs);
        Vector<int>  locs(NProcs * AMREX_SPACEDIM);

        auto comm = ParallelContext::CommunicatorSub();
        ParallelAllGather::AllGather(mx, mxs.dataPtr(), comm);
        BL_ASSERT(sizeof(IntVect) == sizeof(int)*AMREX_SPACEDIM);
        ParallelAllGather::AllGather(loc.getVect(), AMREX_SPACEDIM, locs.dataPtr(), comm);

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

    return loc;
}

Real
MultiFab::norm0 (const iMultiFab& mask, int comp, int nghost, bool local) const
{
    // TODO GPU -- CHECK 

    Real nm0 = -std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(max:nm0) 
#endif
    {
        amrex::Gpu::DeviceScalar<Real> local_mx(-std::numeric_limits<Real>::max());
        Real* p = local_mx.dataPtr();
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            FArrayBox const* fab = &(get(mfi));
            IArrayBox const* mask_fab = &(mask[mfi]);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
            {
                Real t = fab->norminfmask(tbx, *mask_fab, comp, 1);
                amrex::Gpu::Atomic::Max(p, t);
            });
        }
        nm0 = std::max(nm0, local_mx.dataValue());
    }

    if (!local)	ParallelAllReduce::Max(nm0, ParallelContext::CommunicatorSub());

    return nm0;
}

Real
MultiFab::norm0 (int comp, const BoxArray& ba, int nghost, bool local) const
{
    // TODO GPU -- CHECK

    Real nm0 = -std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(max:nm0)
#endif
    {
	std::vector< std::pair<int,Box> > isects;
        amrex::Gpu::DeviceScalar<Real> local_nm0(-std::numeric_limits<Real>::max());
        Real* p = local_nm0.dataPtr();

	for (MFIter mfi(*this); mfi.isValid(); ++mfi)
	{
	    ba.intersections(amrex::grow(mfi.validbox(),nghost),isects);
            FArrayBox const* fab = &(get(mfi));

	    for (int i = 0; i<isects.size(); ++i)
	    {
                const Box& bx = isects[i].second;
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
                {
                    Real t = fab->norm(tbx, 0, comp, 1);
                    amrex::Gpu::Atomic::Max(p, t);
                });

	    }
	}
        nm0 = std::max(nm0, local_nm0.dataValue());
    }
 
    if (!local)
	ParallelAllReduce::Max(nm0, ParallelContext::CommunicatorSub());
 
    return nm0;
}

Real
MultiFab::norm0 (int comp, int nghost, bool local) const
{
    // TODO GPU -- CHECK 

    Real nm0 = -std::numeric_limits<Real>::max();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(max:nm0) 
#endif
    {
        amrex::Gpu::DeviceScalar<Real> local_mx(-std::numeric_limits<Real>::max());
        Real* p = local_mx.dataPtr();
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            FArrayBox const* fab = &(get(mfi));
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
            {
                Real t = fab->norm(tbx, 0, comp, 1);
                amrex::Gpu::Atomic::Max(p, t);
            });
        }
        nm0 = std::max(nm0, local_mx.dataValue());
    }

    if (!local)
	ParallelAllReduce::Max(nm0, ParallelContext::CommunicatorSub());

    return nm0;
}

Vector<Real>
MultiFab::norm0 (const Vector<int>& comps, int nghost, bool local) const
{
    // TODO GPU -- CHECK

    int n = comps.size();
    const Real rmax = std::numeric_limits<Real>::max();
    Vector<Real> nm0(n, -rmax);

#ifdef _OPENMP
    int nthreads = Gpu::notInLaunchRegion() ? omp_get_max_threads() : 1;
#else
    int nthreads = 1;
#endif
    // Device malloc'd copy of comps
    const Gpu::DeviceVector<int> d_comps(comps.begin(), comps.end());
    const int* dptr_comps = d_comps.dataPtr();

    // Host and Device copy of thread's results.
    Vector<Gpu::HostVector<Real> > h_priv_nm0(nthreads, Vector<Real>(n, -rmax));
    Vector<Gpu::DeviceVector<Real> > d_priv_nm0(nthreads, Vector<Real>(n, -rmax));

#ifdef _OPENMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
	int tid = 0;
#endif

        Real* thrd_nm0 = (d_priv_nm0[tid]).dataPtr();
	for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
            const Box& bx = mfi.growntilebox(nghost);
            FArrayBox const* fab = &(get(mfi));
 
            for (int i=0; i<n; i++) {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
                {
                    Real t = fab->norm(tbx, 0, dptr_comps[i], 1);
                    amrex::Gpu::Atomic::Max(&(thrd_nm0[i]), t);
                });
            }
        }

        Gpu::thrust_copy(d_priv_nm0[tid].begin(), d_priv_nm0[tid].end(), h_priv_nm0[tid].begin());

#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
#endif
	for (int i=0; i<n; i++) {
            for (int it=0; it<nthreads; it++) {
	        nm0[i] = std::max(h_priv_nm0[it][i], nm0[i]);
            }	    
	}
    }

    if (!local)
	ParallelAllReduce::Max(nm0.dataPtr(), n, ParallelContext::CommunicatorSub());

    return nm0;
}

Real
MultiFab::norm2 (int comp) const
{
    BL_ASSERT(ixType().cellCentered());

    // Dot expects two MultiFabs. Make a copy to avoid aliasing.
    MultiFab tmpmf(boxArray(), DistributionMap(), 1, 0, MFInfo(), Factory());
    MultiFab::Copy(tmpmf, *this, comp, 0, 1, 0);

    Real nm2 = MultiFab::Dot(*this, comp, tmpmf, 0, 1, 0);
    nm2 = std::sqrt(nm2);
    return nm2;
}

Real
MultiFab::norm2 (int comp, const Periodicity& period) const
{
    MultiFab tmpmf(boxArray(), DistributionMap(), 1, 0, MFInfo(), Factory());
    MultiFab::Copy(tmpmf, *this, comp, 0, 1, 0);

    auto mask = OverlapMask(period);
    MultiFab::Divide(tmpmf, *mask, 0, 0, 1, 0);

    Real nm2 = MultiFab::Dot(*this, comp, tmpmf, 0, 1, 0);
    nm2 = std::sqrt(nm2);
    return nm2;
}

Vector<Real>
MultiFab::norm2 (const Vector<int>& comps) const
{
    // TODO GPU

    BL_ASSERT(ixType().cellCentered());

    int n = comps.size();
    Vector<Real> nm2(n, 0.e0);

#ifdef _OPENMP
    int nthreads = Gpu::notInLaunchRegion() ? omp_get_max_threads() : 1;
#else
    int nthreads = 1;
#endif
    // Device malloc'd copy of comps
    const Gpu::DeviceVector<int> d_comps(comps.begin(), comps.end());
    const int* dptr_comps = d_comps.dataPtr();

    // Host and Device copy of thread's results.
    Vector<Gpu::HostVector<Real> > h_priv_nm2(nthreads, Vector<Real>(n, 0.0));
    Vector<Gpu::DeviceVector<Real> > d_priv_nm2(nthreads, Vector<Real>(n, 0.0));

#ifdef _OPENMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
	int tid = 0;
#endif

        Real* thrd_nm2 = (d_priv_nm2[tid]).dataPtr();
	for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();
	    const FArrayBox* fab = &(get(mfi));
            for (int i=0; i<n; i++) {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
                {
                    Real t = fab->dot(tbx,dptr_comps[i],*fab,tbx,dptr_comps[i]);
                    amrex::Gpu::Atomic::Add((&thrd_nm2[i]), t);
                });
            }
        }

        Gpu::thrust_copy(d_priv_nm2[tid].begin(), d_priv_nm2[tid].end(), h_priv_nm2[tid].begin());

#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
#endif
	for (int i=0; i<n; i++) {
	    for (int it=1; it<nthreads; it++) {
		h_priv_nm2[0][i] += h_priv_nm2[it][i];
	    }
	}
    }

    ParallelAllReduce::Sum(h_priv_nm2[0].dataPtr(), n, ParallelContext::CommunicatorSub());

    for (int i=0; i<n; i++) {
	nm2[i] = std::sqrt(h_priv_nm2[0][i]);
    }

    return nm2;
}

Real
MultiFab::norm1 (int comp, const Periodicity& period) const
{
    MultiFab tmpmf(boxArray(), DistributionMap(), 1, 0, MFInfo(), Factory());
    MultiFab::Copy(tmpmf, *this, comp, 0, 1, 0);

    auto mask = OverlapMask(period);
    MultiFab::Divide(tmpmf, *mask, 0, 0, 1, 0);

    return tmpmf.norm1(0, 0);
}

Real
MultiFab::norm1 (int comp, int ngrow, bool local) const
{
    // TODO GPU -- CHECK
    
    Real nm1 = 0.e0;

#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction && Gpu::notInLaunchRegion()) reduction(+:nm1)
#endif
    {
        amrex::Gpu::DeviceScalar<Real> local_nm1(0.0);
        Real* p = local_nm1.dataPtr();
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(ngrow);
            FArrayBox const* fab = &(get(mfi));
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
            {
                Real t = fab->norm(tbx, 1, comp, 1);
                amrex::Gpu::Atomic::Add(p, t);
            });
        }
        nm1 += local_nm1.dataValue();
    }

    if (!local)
	ParallelAllReduce::Sum(nm1, ParallelContext::CommunicatorSub());

    return nm1;
}

Vector<Real>
MultiFab::norm1 (const Vector<int>& comps, int ngrow, bool local) const
{
    // TODO GPU -- CHECK

    BL_ASSERT(ixType().cellCentered());

    int n = comps.size();
    Vector<Real> nm1(n, 0.e0);

#ifdef _OPENMP
    int nthreads = Gpu::notInLaunchRegion() ? omp_get_max_threads() : 1;
#else
    int nthreads = 1;
#endif
    // Device malloc'd copy of comps
    const Gpu::DeviceVector<int> d_comps(comps.begin(), comps.end());
    const int* dptr_comps = d_comps.dataPtr();

    // Host and Device copy of thread's results.
    Vector<Gpu::HostVector<Real> > h_priv_nm1(nthreads, Vector<Real>(n, 0.0));
    Vector<Gpu::DeviceVector<Real> > d_priv_nm1(nthreads, Vector<Real>(n, 0.0));

#ifdef _OPENMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
	int tid = 0;
#endif

        Real* thrd_nm1 = (d_priv_nm1[tid]).dataPtr();
	for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
            const Box& bx = mfi.growntilebox(ngrow);
            FArrayBox const* fab = &(get(mfi));

            for (int i=0; i<n; i++) {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
                {
                    Real t = fab->norm(tbx, 1, dptr_comps[i], 1);
                    amrex::Gpu::Atomic::Max(&(thrd_nm1[i]), t);
                });
	    }
        }

        Gpu::thrust_copy(d_priv_nm1[tid].begin(), d_priv_nm1[tid].end(), h_priv_nm1[tid].begin());

#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
#endif
	for (int i=0; i<n; i++) {
	    for (int it=0; it<nthreads; it++) {
		nm1[i] += h_priv_nm1[it][i];
	    }
	}
    }

    if (!local)
	ParallelAllReduce::Sum(nm1.dataPtr(), n, ParallelContext::CommunicatorSub());

    return nm1;
}

Real
MultiFab::sum (int comp, bool local) const
{
    // TODO GPU -- CHECK

    Real sm = 0.e0;

#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction && Gpu::notInLaunchRegion()) reduction(+:sm)
#endif
    {
        amrex::Gpu::DeviceScalar<Real> local_sm(0.0);
        Real* p = local_sm.dataPtr();
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            FArrayBox const* fab = &(get(mfi));
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
            {
                Real t = fab->sum(tbx, comp, 1);
                amrex::Gpu::Atomic::Add(p, t);
            });
        }
        sm += local_sm.dataValue();
    }

    if (!local)
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());

    return sm;
}

void
MultiFab::minus (const MultiFab& mf,
                 int             strt_comp,
                 int             num_comp,
                 int             nghost)
{
    MultiFab::Subtract(*this, mf, strt_comp, strt_comp, num_comp, nghost);
}

void
MultiFab::divide (const MultiFab& mf,
		  int             strt_comp,
		  int             num_comp,
		  int             nghost)
{
    MultiFab::Divide(*this, mf, strt_comp, strt_comp, num_comp, nghost);
}

void
MultiFab::plus (Real val,
                int  comp,
                int  num_comp,
                int  nghost)
{

    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        FArrayBox* fab = &(get(mfi));
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
        {
            fab->plus(val, tbx, comp, num_comp);
        });
    }
}

void
MultiFab::plus (Real       val,
                const Box& region,
                int        comp,
                int        num_comp,
                int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& b = mfi.growntilebox(nghost) & region;
        FArrayBox* fab = &(get(mfi));
        if (b.ok()) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( b, tbx,
            {
                fab->plus(val, tbx, comp, num_comp);
            });
        }
    }
}

void
MultiFab::plus (const MultiFab& mf,
                int             strt_comp,
                int             num_comp,
                int             nghost)
{
    MultiFab::Add(*this, mf, strt_comp, strt_comp, num_comp, nghost);
}

void
MultiFab::mult (Real val,
                int  comp,
                int  num_comp,
                int  nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        FArrayBox* fab = &(get(mfi));
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
        {
            fab->mult(val, tbx, comp, num_comp);
        });
    }
}

void
MultiFab::mult (Real       val,
                const Box& region,
                int        comp,
                int        num_comp,
                int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& b = mfi.growntilebox(nghost) & region;

        if (b.ok()) {
            FArrayBox* fab = &(get(mfi));
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( b, tbx,
            {
                fab->mult(val, tbx, comp, num_comp);
            });
        }
    }
}

void
MultiFab::invert (Real numerator,
                  int  comp,
                  int  num_comp,
                  int  nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        FArrayBox* fab = &(get(mfi));
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
        {
            fab->invert(numerator, tbx, comp, num_comp);
        });
    }
}

void
MultiFab::invert (Real       numerator,
                  const Box& region,
                  int        comp,
                  int        num_comp,
                  int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);
    BL_ASSERT(num_comp > 0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& b = mfi.growntilebox(nghost) & region;

        if (b.ok()) {
            FArrayBox* fab = &(get(mfi));
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( b, tbx,
            {
                fab->invert(numerator, tbx, comp, num_comp);
            });
        }
    }
}

void
MultiFab::negate (int comp,
                  int num_comp,
                  int nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        FArrayBox* fab = &(get(mfi));
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
        {
            fab->negate(tbx, comp, num_comp);
        });
    }
}

void
MultiFab::negate (const Box& region,
                  int        comp,
                  int        num_comp,
                  int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow.min());
    BL_ASSERT(comp+num_comp <= n_comp);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& b = mfi.growntilebox(nghost) & region;

        if (b.ok()) {
            FArrayBox* fab = &(get(mfi));
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA( b, tbx,
            {
                fab->negate(tbx, comp, num_comp);
            });
        }
    }
}

void
MultiFab::SumBoundary (int scomp, int ncomp, const Periodicity& period)
{
    BL_PROFILE("MultiFab::SumBoundary()");

    if ( n_grow == IntVect::TheZeroVector() && boxArray().ixType().cellCentered()) return;

    if (boxArray().ixType().cellCentered()) {
	// Self copy is safe only for cell-centered MultiFab
	this->copy(*this,scomp,scomp,ncomp,n_grow,IntVect::TheZeroVector(),period,FabArrayBase::ADD);
    } else {
	MultiFab tmp(boxArray(), DistributionMap(), ncomp, n_grow, MFInfo(), Factory());
	MultiFab::Copy(tmp, *this, scomp, 0, ncomp, n_grow);
	this->setVal(0.0, scomp, ncomp, 0);
	this->copy(tmp,0,scomp,ncomp,n_grow,IntVect::TheZeroVector(),period,FabArrayBase::ADD);
    }
}

void
MultiFab::SumBoundary (const Periodicity& period)
{
    SumBoundary(0, n_comp, period);
}

std::unique_ptr<MultiFab>
MultiFab::OverlapMask (const Periodicity& period) const
{
    //TODO GPU????

    BL_PROFILE("MultiFab::OverlapMask()");

    const BoxArray& ba = boxArray();
    const DistributionMapping& dm = DistributionMap();

    std::unique_ptr<MultiFab> p{new MultiFab(ba,dm,1,0, MFInfo(), Factory())};
    p->setVal(0.0);

    const std::vector<IntVect>& pshifts = period.shiftIntVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector< std::pair<int,Box> > isects;
        
        for (MFIter mfi(*p); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = (*p)[mfi];
            const Box& bx = fab.box();
            for (const auto& iv : pshifts)
            {
                ba.intersections(bx+iv, isects);                    
                for (const auto& is : isects)
                {
                    fab.plus(1.0, is.second-iv);
                }
            }
        }
    }
    
    return p;
}

std::unique_ptr<iMultiFab>
MultiFab::OwnerMask (const Periodicity& period) const
{
    //TODO GPU????
    BL_PROFILE("MultiFab::OwnerMask()");

    const BoxArray& ba = boxArray();
    const DistributionMapping& dm = DistributionMap();

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
    tmpmf.setVal(0.0);
    tmpmf.ParallelCopy(*this, period, FabArrayBase::ADD);

    MultiFab::Copy(*this, tmpmf, 0, 0, ncomp, 0);
}

void
MultiFab::OverrideSync (const Periodicity& period)
{
    if (ixType().cellCentered()) return;
    auto msk = this->OwnerMask(period);
    this->OverrideSync(*msk, period);
}

void
MultiFab::OverrideSync (const iMultiFab& msk, const Periodicity& period)
{
    BL_PROFILE("MultiFab::OverrideSync()");

    if (ixType().cellCentered()) return;
    
    const int ncomp = nComp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox* fab = &(*this)[mfi];
        IArrayBox const* ifab = &(msk[mfi]);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA( bx, tbx,
        {
            fab->setValIfNot(0.0, tbx, *ifab, 0, ncomp);
        });
    }
    
    MultiFab tmpmf(boxArray(), DistributionMap(), ncomp, 0, MFInfo(), Factory());
    tmpmf.setVal(0.0);
    tmpmf.ParallelCopy(*this, period, FabArrayBase::ADD);

    MultiFab::Copy(*this, tmpmf, 0, 0, ncomp, 0);
}

}
