
#include <FabSet.H>
#include <ParallelDescriptor.H>
#include <BLProfiler.H>
#include <VisMF.H>

#ifdef _OPENMP
#include <omp.h>
#endif

FabSet::FabSet () {}

FabSet::~FabSet () {}

FabSet::FabSet (const BoxArray& grids, int ncomp, ParallelDescriptor::Color color)
    :
    m_mf(grids,ncomp,0,color)
{}

void
FabSet::define (const BoxArray& grids, int ncomp, ParallelDescriptor::Color color)
{
    m_mf.define(grids, ncomp, 0, Fab_allocate, IntVect::TheZeroVector(), color);
}

void
FabSet::define (const BoxArray& grids, int ncomp, const DistributionMapping& dm)
{
    m_mf.define(grids, ncomp, 0, dm, Fab_allocate);
}

FabSet&
FabSet::copyFrom (const FabSet& src, int scomp, int dcomp, int ncomp)
{
    if (boxArray() == src.boxArray() && DistributionMap() == src.DistributionMap()) {
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (FabSetIter fsi(*this); fsi.isValid(); ++fsi) {
	    (*this)[fsi].copy(src[fsi], scomp, dcomp, ncomp);
	}
    } else {
	m_mf.copy(src.m_mf,scomp,dcomp,ncomp);
    }
    return *this;
}

FabSet&
FabSet::copyFrom (const MultiFab& src, int ngrow, int scomp, int dcomp, int ncomp)
{
    BL_ASSERT(boxArray() != src.boxArray());
    m_mf.copy(src,scomp,dcomp,ncomp,ngrow,0);
    return *this;
}

FabSet&
FabSet::plusFrom (const FabSet& src, int scomp, int dcomp, int ncomp)
{
    if (boxArray() == src.boxArray() && DistributionMap() == src.DistributionMap()) {
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (FabSetIter fsi(*this); fsi.isValid(); ++fsi) {
	    (*this)[fsi].plus(src[fsi], scomp, dcomp, ncomp);
	}
    } else {
	BoxLib::Abort("FabSet::plusFrom: parallel plusFrom not supported");
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
FabSet::setVal (Real val)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi) {
	(this->m_mf)[fsi].setVal(val);
    }
}

void
FabSet::setVal (Real val, int comp, int num_comp)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi) {
	FArrayBox& fab = (this->m_mf)[fsi];
	fab.setVal(val, fab.box(), comp, num_comp);
    }
}

// Linear combination this := a*this + b*src
// Note: corresponding fabsets must be commensurate.
FabSet&
FabSet::linComb (Real a, Real b, const FabSet& src, int scomp, int dcomp, int ncomp)
{
    BL_ASSERT(size() == src.size());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
	const FArrayBox& srcfab = src[fsi];
	FArrayBox& dstfab = (*this)[fsi];
	BL_ASSERT(srcfab.box() == dstfab.box());
	dstfab.mult(a, dcomp, ncomp);
	dstfab.saxpy(b, srcfab, srcfab.box(), dstfab.box(), scomp, dcomp, ncomp);
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

    MultiFab bdrya(boxArray(),ncomp,0,DistributionMap());
    MultiFab bdryb(boxArray(),ncomp,0,DistributionMap());

    bdrya.copy(mfa,a_comp,0,ncomp,ngrow,0);
    bdryb.copy(mfb,b_comp,0,ncomp,ngrow,0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
	const FArrayBox& afab = bdrya[fsi];
	const FArrayBox& bfab = bdryb[fsi];
	FArrayBox& dfab = (*this)[fsi];
	dfab.linComb(afab, afab.box(), a_comp,
		     bfab, bfab.box(), b_comp,
		     a, b, dfab.box(), dcomp, ncomp);
    }

    return *this;
}

void
FabSet::write(const std::string& name) const
{
    // xxxxx FIXME
    VisMF::Write(m_mf,name);
}

void
FabSet::read(const std::string& name)
{
    //xxxxx FIXME
    VisMF::Read(m_mf,name);
}
