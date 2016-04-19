
#include <FillPatchUtil.H>

#include <iostream>

namespace BoxLib
{
    void FillPatchSingleLevel (MultiFab& mf, Real time, const PArray<MultiFab>& smf,
			       const std::vector<Real>& stime, int scomp,
			       const Geometry& geom, PhysBCFunct& physbcf)
    {
	int ncomp = mf.nComp();

	BL_ASSERT(scomp+ncomp <= smf[0].nComp());
	BL_ASSERT(smf.size() == stime.size());

	if (smf.size() == 1) 
	{
	    mf.copy(smf[0], scomp, 0, ncomp);

	    mf.FillBoundary_nowait();
	    geom.FillPeriodicBoundary_nowait(mf);
	
	    mf.FillBoundary_finish();
	    geom.FillPeriodicBoundary_finish(mf);

	    physbcf.doit(mf, time);
	} 
	else if (smf.size() == 2) 
	{
	    std::cout << "xxxxx " << time << ", " << stime[0] << ", " << stime[1] << std::endl;

	    BoxLib::Abort("FillPatchSingleLevel: linear interpolation in time not implemented yet");
	} else {
	    BoxLib::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
	}
    }
}
