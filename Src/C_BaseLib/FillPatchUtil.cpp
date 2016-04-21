
#include <FillPatchUtil.H>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace BoxLib
{
    void FillPatchSingleLevel (MultiFab& mf, Real time, 
			       const PArray<MultiFab>& smf, const std::vector<Real>& stime,
			       int scomp, int dcomp, int ncomp,
			       const Geometry& geom, PhysBCFunct& physbcf)
    {
	BL_PROFILE("FillPatchSingleLevel");

	BL_ASSERT(scomp+ncomp <= smf[0].nComp());
	BL_ASSERT(dcomp+ncomp <= mf.nComp());
	BL_ASSERT(smf.size() == stime.size());
	BL_ASSERT(smf.size() != 0);

	if (smf.size() <= 2) 
	{
	    if (smf.size() == 1) {
		mf.copy(smf[0], scomp, dcomp, ncomp);
		geom.PeriodicCopy(mf, smf[0], dcomp, scomp, ncomp);
	    } else {
		BL_ASSERT(smf[0].boxArray() == smf[1].boxArray());
		PArray<MultiFab> raii(PArrayManage);
		MultiFab * dmf;
		int destcomp;
		bool sameba;
		if (mf.boxArray() == smf[0].boxArray()) {
		    dmf = &mf;
		    destcomp = dcomp;
		    sameba = true;
		} else {
		    dmf = raii.push_back(new MultiFab(smf[0].boxArray(), ncomp, 0));
		    destcomp = 0;
		    sameba = false;
		}

#ifdef _OPENMP
#pragma omp parallel 
#endif
		for (MFIter mfi(*dmf,true); mfi.isValid(); ++mfi)
		{
		    const Box& bx = mfi.tilebox();
		    (*dmf)[mfi].linInterp(smf[0][mfi],
					  scomp,
					  smf[1][mfi],
					  scomp,
					  stime[0],
					  stime[1],
					  time,
					  bx,
					  destcomp,
					  ncomp);
		}

		if (!sameba) {
		    int src_ngrow = 0;
		    int dst_ngrow = mf.nGrow();
		    mf.copy(*dmf, 0, dcomp, ncomp, src_ngrow, dst_ngrow);
		    geom.PeriodicCopy(mf, *dmf, dcomp, 0, ncomp);
		}
	    }

	    mf.FillBoundary_nowait(dcomp,ncomp);
	    geom.FillPeriodicBoundary_nowait(mf,dcomp,ncomp);
	
	    mf.FillBoundary_finish();
	    geom.FillPeriodicBoundary_finish(mf);

	    physbcf.doit(mf, dcomp, ncomp, time);
	} 
	else {
	    BoxLib::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
	}
    }


    void FillPatchTwoLevels (MultiFab& mf, Real time,
			     const PArray<MultiFab>& cmf, const std::vector<Real>& ct,
			     const PArray<MultiFab>& fmf, const std::vector<Real>& ft,
			     int scomp, int dcomp, int ncomp,
			     const Geometry& cgeom, const Geometry& fgeom, 
			     PhysBCFunct& cbc, PhysBCFunct& fbc,
			     const IntVect& ratio, 
			     Interpolater* mapper, const Array<BCRec>& bcs)
    {
	const BoxArray&  ba =     mf.boxArray();
	const BoxArray& fba = fmf[0].boxArray();

	const DistributionMapping& dm = mf.DistributionMap();

	const IndexType& boxtype = ba.ixType();

	const int myproc = ParallelDescriptor::MyProc();

	const int ngrow = mf.nGrow();

	Box fdomain = fgeom.Domain();
	fdomain.convert(boxtype);
	Box fdomain_g(fdomain);
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    if (fgeom.isPeriodic(i)) {
		fdomain_g.grow(i,ngrow);
	    }
	}

	BoxList bl(boxtype);
	Array<int> idxs;
	Array<int> iprocs;
	Array<Box> fpatch;
	std::vector< std::pair<int,Box> > isects;

	for (int i = 0, N = ba.size(); i < N; ++i)
	{
	    Box bx = ba[i];
	    bx.grow(ngrow);
	    bx &= fdomain_g;

	    fba.intersections(bx, isects);

	    BoxList pieces(boxtype);
	    for (std::vector< std::pair<int,Box> >::const_iterator it = isects.begin();
		 it != isects.end(); ++it)
	    {
		pieces.push_back(it->second);
	    }
	    BoxList leftover = BoxLib::complementIn(bx, pieces);

	    bool ismybox = (dm[i] == myproc);
	    for (BoxList::const_iterator bli = leftover.begin(); bli != leftover.end(); ++bli)
	    {
		bl.push_back(mapper->CoarseBox(*bli,ratio));
		if (ismybox) {
		    fpatch.push_back(*bli);
		    idxs.push_back(i);
		}
		iprocs.push_back(dm[i]);
	    }
	}
	
	if (!iprocs.empty())
	{
	    BoxArray ba_crse_patch(bl);

	    iprocs.push_back(ParallelDescriptor::MyProc());
	    DistributionMapping dm_crse_patch(iprocs,false);

	    MultiFab mf_crse_patch(ba_crse_patch, ncomp, 0, dm_crse_patch);

	    FillPatchSingleLevel(mf_crse_patch, time, cmf, ct, scomp, 0, ncomp, cgeom, cbc);

	    int idummy1=0, idummy2=0;
#ifdef _OPENMP
#pragma omp parallel
#endif
	    for (MFIter mfi(mf_crse_patch); mfi.isValid(); ++mfi)
	    {
		int li = mfi.LocalIndex();
		int gi = idxs[li];		
		const Box& dbx = fpatch[li];

		Array<BCRec> bcr(ncomp);
		BoxLib::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);
		
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
			       idummy1, idummy2);
	    }
	}

	FillPatchSingleLevel(mf, time, fmf, ft, scomp, dcomp, ncomp, fgeom, fbc);
    }
}
