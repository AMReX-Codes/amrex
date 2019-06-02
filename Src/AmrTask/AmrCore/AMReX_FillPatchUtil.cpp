#include <AMReX_Utility.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FillPatchUtil_F.H>
#include <cmath>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex
{
    void FillPatchSingleLevel (MultiFab& mf, Real time, 
	    const Vector<MultiFab*>& smf, const Vector<Real>& stime,
	    int scomp, int dcomp, int ncomp,
	    const Geometry& geom, PhysBCFunctBase& physbcf)
    {
	BL_PROFILE("FillPatchSingleLevel");

	BL_ASSERT(scomp+ncomp <= smf[0]->nComp());
	BL_ASSERT(dcomp+ncomp <= mf.nComp());
	BL_ASSERT(smf.size() == stime.size());
	BL_ASSERT(smf.size() != 0);

	if (smf.size() == 1) 
	{
	    mf.copy(*smf[0], scomp, dcomp, ncomp, 0, mf.nGrow(), geom.periodicity());
	} 
	else if (smf.size() == 2) 
	{
	    BL_ASSERT(smf[0]->boxArray() == smf[1]->boxArray());
	    MultiFab raii;
	    MultiFab * dmf;
	    int destcomp;
	    bool sameba;
	    if (mf.boxArray() == smf[0]->boxArray()) {
		dmf = &mf;
		destcomp = dcomp;
		sameba = true;
	    } else {
		raii.define(smf[0]->boxArray(), smf[0]->DistributionMap(), ncomp, 0,
			MFInfo(), smf[0]->Factory());

		dmf = &raii;
		destcomp = 0;
		sameba = false;
	    }

#ifdef _OPENMP
#pragma omp parallel 
#endif
	    for (MFIter mfi(*dmf,true); mfi.isValid(); ++mfi)
	    {
		const Box& bx = mfi.tilebox();
		(*dmf)[mfi].linInterp((*smf[0])[mfi],
			scomp,
			(*smf[1])[mfi],
			scomp,
			stime[0],
			stime[1],
			time,
			bx,
			destcomp,
			ncomp);
	    }

	    if (sameba)
	    {
		// Note that when sameba is true mf's BoxArray is nonoverlapping.
		// So FillBoundary is safe.
		mf.FillBoundary(dcomp,ncomp,geom.periodicity());
	    }
	    else
	    {
		int src_ngrow = 0;
		int dst_ngrow = mf.nGrow();

		mf.copy(*dmf, 0, dcomp, ncomp, src_ngrow, dst_ngrow, geom.periodicity());
	    }
	}
	else {
	    amrex::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
	}

	physbcf.FillBoundary(mf, dcomp, ncomp, time);
    }

    void FillPatchTwoLevels (MultiFab& mf, Real time,
	    const Vector<MultiFab*>& cmf, const Vector<Real>& ct,
	    const Vector<MultiFab*>& fmf, const Vector<Real>& ft,
	    int scomp, int dcomp, int ncomp,
	    const Geometry& cgeom, const Geometry& fgeom, 
	    PhysBCFunctBase& cbc, PhysBCFunctBase& fbc,
	    const IntVect& ratio, 
	    Interpolater* mapper, const BCRec& bcs)
    {
	Vector<BCRec> bcs_array(1,BCRec(bcs.lo(),bcs.hi()));

	FillPatchTwoLevels(mf,time,cmf,ct,fmf,ft,scomp,dcomp,ncomp,cgeom,fgeom,
		cbc,fbc,ratio,mapper,bcs_array);
    }


    void FillPatchTwoLevels (MultiFab& mf, Real time,
	    const Vector<MultiFab*>& cmf, const Vector<Real>& ct,
	    const Vector<MultiFab*>& fmf, const Vector<Real>& ft,
	    int scomp, int dcomp, int ncomp,
	    const Geometry& cgeom, const Geometry& fgeom, 
	    PhysBCFunctBase& cbc, PhysBCFunctBase& fbc,
	    const IntVect& ratio, 
	    Interpolater* mapper, const Vector<BCRec>& bcs)
    {
	BL_PROFILE("FillPatchTwoLevels");

	int ngrow = mf.nGrow();

	if (ngrow > 0 || mf.getBDKey() != fmf[0]->getBDKey()) 
	{
	    const InterpolaterBoxCoarsener& coarsener = mapper->BoxCoarsener(ratio);

	    Box fdomain = fgeom.Domain();
	    fdomain.convert(mf.boxArray().ixType());
	    Box fdomain_g(fdomain);
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
		if (fgeom.isPeriodic(i)) {
		    fdomain_g.grow(i,ngrow);
		}
	    }

	    const FabArrayBase::FPinfo& fpc = FabArrayBase::TheFPinfo(*fmf[0], mf, fdomain_g,
                                                                      IntVect(ngrow), coarsener, 
                                                                      amrex::coarsen(fgeom.Domain(),ratio));

	    if ( ! fpc.ba_crse_patch.empty())
	    {
		MultiFab mf_crse_patch(fpc.ba_crse_patch, fpc.dm_crse_patch, ncomp, 0, MFInfo(),
			*fpc.fact_crse_patch);

		FillPatchSingleLevel(mf_crse_patch, time, cmf, ct, scomp, 0, ncomp, cgeom, cbc);

		int idummy1=0, idummy2=0;
		bool cc = fpc.ba_crse_patch.ixType().cellCentered();
		ignore_unused(cc);
#ifdef _OPENMP
#pragma omp parallel if (cc)
#endif
		for (MFIter mfi(mf_crse_patch); mfi.isValid(); ++mfi)
		{
		    int li = mfi.LocalIndex();
		    int gi = fpc.dst_idxs[li];		
		    const Box& dbx = fpc.dst_boxes[li];

		    Vector<BCRec> bcr(ncomp);
		    amrex::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);

		    mapper->interp(mf_crse_patch[mfi],
			    0,
			    mf[gi],
			    dcomp,
			    ncomp,
			    dbx,
			    ratio,
			    cgeom,
			    fgeom,
			    bcr,
			    idummy1, idummy2, RunOn::Cpu);
		}
	    }
	}

	FillPatchSingleLevel(mf, time, fmf, ft, scomp, dcomp, ncomp, fgeom, fbc);
    }
}//namespace
