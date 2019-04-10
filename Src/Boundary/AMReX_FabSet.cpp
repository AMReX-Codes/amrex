
#include <AMReX_FabSet.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_VisMF.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

FabSet::FabSet () noexcept {}

FabSet::FabSet (const BoxArray& grids, const DistributionMapping& dmap, int ncomp)
    :
    m_mf(grids,dmap,ncomp,0,MFInfo(),FArrayBoxFactory())
{}

void
FabSet::define (const BoxArray& grids, const DistributionMapping& dm, int ncomp)
{
    m_mf.define(grids, dm, ncomp, 0, MFInfo(), FArrayBoxFactory());
}

FabSet&
FabSet::copyFrom (const FabSet& src, int scomp, int dcomp, int ncomp)
{
    if (boxArray() == src.boxArray() && DistributionMap() == src.DistributionMap()) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for (FabSetIter fsi(*this); fsi.isValid(); ++fsi) {
            const Box& bx = fsi.validbox();
            auto const srcfab =   src.array(fsi);
            auto       dstfab = this->array(fsi);
            AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                dstfab(i,j,k,n+dcomp) = srcfab(i,j,k,n+scomp);
            });
	}
    } else {
	m_mf.copy(src.m_mf,scomp,dcomp,ncomp);
    }
    return *this;
}

FabSet&
FabSet::copyFrom (const MultiFab& src, int ngrow, int scomp, int dcomp, int ncomp,
		  const Periodicity& period)
{
    BL_ASSERT(boxArray() != src.boxArray());
    m_mf.copy(src,scomp,dcomp,ncomp,ngrow,0,period);
    return *this;
}

FabSet&
FabSet::plusFrom (const FabSet& src, int scomp, int dcomp, int ncomp)
{
    if (boxArray() == src.boxArray() && DistributionMap() == src.DistributionMap()) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for (FabSetIter fsi(*this); fsi.isValid(); ++fsi) {
            const Box& bx = fsi.validbox();
            auto const srcfab =   src.array(fsi);
            auto       dstfab = this->array(fsi);
            AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                dstfab(i,j,k,n+dcomp) += srcfab(i,j,k,n+scomp);
            });
	}
    } else {
	amrex::Abort("FabSet::plusFrom: parallel plusFrom not supported");
    }
    return *this;
}

FabSet&
FabSet::plusFrom (const MultiFab& src, int ngrow, int scomp, int dcomp, int ncomp,
		  const Periodicity& period)
{
    BL_ASSERT(boxArray() != src.boxArray());
    m_mf.copy(src,scomp,dcomp,ncomp,ngrow,0,period,FabArrayBase::ADD);
    return *this;
}

void
FabSet::copyTo (MultiFab& dest, int ngrow, int scomp, int dcomp, int ncomp,
		const Periodicity& period) const
{
    BL_ASSERT(boxArray() != dest.boxArray());
    dest.copy(m_mf,scomp,dcomp,ncomp,0,ngrow,period);
}

void
FabSet::plusTo (MultiFab& dest, int ngrow, int scomp, int dcomp, int ncomp,
		const Periodicity& period) const
{
    BL_ASSERT(boxArray() != dest.boxArray());
    dest.copy(m_mf,scomp,dcomp,ncomp,0,ngrow,period,FabArrayBase::ADD);
}

void
FabSet::setVal (Real val)
{
    const int ncomp = nComp();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi) {
        const Box& bx = fsi.validbox();
        auto fab = this->array(fsi);
        AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            fab(i,j,k,n) = val;
        });
    }
}

void
FabSet::setVal (Real val, int comp, int num_comp)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi) {
        const Box& bx = fsi.validbox();
        auto fab = this->array(fsi);
        AMREX_HOST_DEVICE_FOR_4D ( bx, num_comp, i, j, k, n,
        {
            fab(i,j,k,n+comp) = val;
        });
    }
}

// Linear combination this := a*this + b*src
// Note: corresponding fabsets must be commensurate.
FabSet&
FabSet::linComb (Real a, Real b, const FabSet& src, int scomp, int dcomp, int ncomp)
{
    BL_ASSERT(size() == src.size());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
        const Box& bx = fsi.validbox();
        auto const srcfab =   src.array(fsi);
        auto       dstfab = this->array(fsi);
        AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            dstfab(i,j,k,n+dcomp) = a*dstfab(i,j,k,n+dcomp) + b*srcfab(i,j,k,n+scomp);
        });
    }
    return *this;
}

// Linear combination: this := a*mfa + b*mfb
// CastroRadiation is the only code that uses this function. 
FabSet&
FabSet::linComb (Real a, const MultiFab& mfa, int a_comp,
		 Real b, const MultiFab& mfb, int b_comp,
		 int dcomp, int ncomp, int ngrow)
{
    BL_PROFILE("FabSet::linComb()");
    BL_ASSERT(ngrow <= mfa.nGrow());
    BL_ASSERT(ngrow <= mfb.nGrow());
    BL_ASSERT(mfa.boxArray() == mfb.boxArray());
    BL_ASSERT(boxArray() != mfa.boxArray());

    MultiFab bdrya(boxArray(),DistributionMap(),ncomp,0,MFInfo(),FArrayBoxFactory());
    MultiFab bdryb(boxArray(),DistributionMap(),ncomp,0,MFInfo(),FArrayBoxFactory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(bdrya); mfi.isValid(); ++mfi) // tiling is not safe for this BoxArray
    {
        const Box& bx = mfi.validbox();
        auto afab = bdrya.array(mfi);
        auto bfab = bdryb.array(mfi);
        AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            afab(i,j,k,n) = 1.e200;
            bfab(i,j,k,n) = 1.e200;
        });
    }

    bdrya.copy(mfa,a_comp,0,ncomp,ngrow,0);
    bdryb.copy(mfb,b_comp,0,ncomp,ngrow,0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
        const Box& bx = fsi.validbox();
        auto const afab = bdrya.array(fsi);
        auto const bfab = bdryb.array(fsi);
        auto       dfab = this->array(fsi);
        AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            dfab(i,j,k,n+dcomp) = a*afab(i,j,k,n) + b*bfab(i,j,k,n);
        });
    }

    return *this;
}

void
FabSet::write(const std::string& name) const
{
    VisMF::Write(m_mf,name);
}

void
FabSet::read(const std::string& name)
{
    if (m_mf.empty()) {
	amrex::Abort("FabSet::read: not predefined");
    }
    VisMF::Read(m_mf,name);
}

void
FabSet::Copy (FabSet& dst, const FabSet& src)
{
    BL_ASSERT(amrex::match(dst.boxArray(), src.boxArray()));
    BL_ASSERT(dst.DistributionMap() == src.DistributionMap());
    int ncomp = dst.nComp();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (FabSetIter fsi(dst); fsi.isValid(); ++fsi) {
        const Box& bx = fsi.validbox();
        auto const srcfab = src.array(fsi);
        auto       dstfab = dst.array(fsi);
        AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            dstfab(i,j,k,n) = srcfab(i,j,k,n);
        });
    }
}

}
