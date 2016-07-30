
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
    m_mf.copy(src,scomp,dcomp,ncomp,ngrow,0,FabArrayBase::COPY);
    return *this;
}

FabSet&
FabSet::plusFrom (const MultiFab& src, int ngrow, int scomp, int dcomp, int ncomp)
{
    BL_ASSERT(boxArray() != src.boxArray());
    m_mf.copy(src,scomp,dcomp,ncomp,ngrow,0,FabArrayBase::ADD);
    return *this;
}

void
FabSet::copyTo (MultiFab& dest) const
{
    BL_ASSERT(boxArray() != dest.boxArray());
    dest.copy(m_mf);
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
        //
        // WARNING: same fab used as src and dest here.
	// xxxxx FIXME
        //
	const FArrayBox& srcfab = src[fsi];
	FArrayBox& dstfab = (*this)[fsi];
	BL_ASSERT(srcfab.box() == dstfab.box());
	dstfab.linComb(dstfab, dstfab.box(), dcomp,
		       srcfab, srcfab.box(), scomp,
		       a, b, dstfab.box(), dcomp, ncomp);
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

#if 0

    BL_ASSERT(ngrow <= mfa.nGrow());
    BL_ASSERT(ngrow <= mfb.nGrow());

    const BoxArray& bxa = mfa.boxArray();

    BL_ASSERT(bxa == mfb.boxArray());

    MultiFabCopyDescriptor mfcd;

    MultiFabId mfid_mfa = mfcd.RegisterFabArray(const_cast<MultiFab*>(&mfa));
    MultiFabId mfid_mfb = mfcd.RegisterFabArray(const_cast<MultiFab*>(&mfb));

    std::vector<FillBoxId> fbids_mfa, fbids_mfb;

    std::vector< std::pair<int,Box> > isects;

    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
        bxa.intersections(get(fsi).box(),isects,ngrow);

        const int index = fsi.index();

        for (int j = 0, N = isects.size(); j < N; j++)
        {
            const int  grd  = isects[j].first;
            const Box& ovlp = isects[j].second;

            fbids_mfa.push_back(mfcd.AddBox(mfid_mfa,
                                            ovlp,
                                            0,
                                            grd,
                                            a_comp,
                                            0,
                                            ncomp,
                                            false));

            BL_ASSERT(fbids_mfa.back().box() == ovlp);
            //
            // Also save the index of the FAB in the FabSet.
            //
            fbids_mfa.back().FabIndex(index);

            fbids_mfb.push_back(mfcd.AddBox(mfid_mfb,
                                            ovlp,
                                            0,
                                            grd,
                                            b_comp,
                                            0,
                                            ncomp,
                                            false));

            BL_ASSERT(fbids_mfb.back().box() == ovlp);
        }
    }

    BL_COMM_PROFILE_NAMETAG("CD::FabSet::linComb()");
    mfcd.CollectData();

    FArrayBox a_fab, b_fab;

    BL_ASSERT(fbids_mfa.size() == fbids_mfb.size());

    for (int i = 0, N = fbids_mfa.size(); i < N; i++)
    {
        a_fab.resize(fbids_mfa[i].box(), ncomp);
        b_fab.resize(fbids_mfb[i].box(), ncomp);

        mfcd.FillFab(mfid_mfa, fbids_mfa[i], a_fab);
        mfcd.FillFab(mfid_mfb, fbids_mfb[i], b_fab);

        BL_ASSERT(DistributionMap()[fbids_mfa[i].FabIndex()] == ParallelDescriptor::MyProc());

        (*this)[fbids_mfa[i].FabIndex()].linComb(a_fab,
                                                 fbids_mfa[i].box(),
                                                 0,
                                                 b_fab,
                                                 fbids_mfa[i].box(),
                                                 0,
                                                 a,
                                                 b,
                                                 fbids_mfa[i].box(),
                                                 dcomp,
                                                 ncomp);
    }
#endif

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
