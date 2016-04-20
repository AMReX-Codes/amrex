
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
}
