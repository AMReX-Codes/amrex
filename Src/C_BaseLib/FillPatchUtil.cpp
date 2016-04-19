
#include <FillPatchUtil.H>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace BoxLib
{
    void FillPatchSingleLevel (MultiFab& mf, Real time, const PArray<MultiFab>& smf,
			       const std::vector<Real>& stime, int scomp,
			       const Geometry& geom, PhysBCFunct& physbcf)
    {
	BL_PROFILE("FillPatchSingleLevel");

	int ncomp = mf.nComp();

	BL_ASSERT(scomp+ncomp <= smf[0].nComp());
	BL_ASSERT(smf.size() == stime.size());

	if (smf.size() <= 2) 
	{
	    if (smf.size() == 1) {
		mf.copy(smf[0], scomp, 0, ncomp);
	    } else {
		BL_ASSERT(smf[0].boxArray() == smf[1].boxArray());
		if (mf.boxArray() == smf[0].boxArray()) {
#ifdef _OPENMP
#pragma omp parallel 
#endif
		    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
		    {
			const Box& bx = mfi.tilebox();
			mf[mfi].linInterp(smf[0][mfi],
					  scomp,
					  smf[1][mfi],
					  scomp,
					  stime[0],
					  stime[1],
					  time,
					  bx,
					  0,
					  ncomp);
		    }
		} else {
		    BoxLib::Abort("FillPatchSingleLevel: MFs do not have the same BoxArray");
		}
	    }

	    mf.FillBoundary_nowait();
	    geom.FillPeriodicBoundary_nowait(mf);
	
	    mf.FillBoundary_finish();
	    geom.FillPeriodicBoundary_finish(mf);

	    physbcf.doit(mf, time);
	} 
	else {
	    BoxLib::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
	}
    }
}
