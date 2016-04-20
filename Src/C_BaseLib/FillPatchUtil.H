#ifndef _FillPatchUtil_H_
#define _FillPatchUtil_H_

#include <MultiFab.H>
#include <Geometry.H>
#include <PhysBCFunct.H>

namespace BoxLib
{
    void FillPatchSingleLevel (MultiFab& mf, Real time, 
			       const PArray<MultiFab>& smf, const std::vector<Real>& stime, 
			       int scomp, int dcomp, int ncomp,
			       const Geometry& geom, PhysBCFunct& physbcf);
}

#endif
