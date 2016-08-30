#ifndef _FillPatchUtil_H_
#define _FillPatchUtil_H_

#include <MultiFab.H>
#include <Geometry.H>
#include <PhysBCFunct.H>
#include <Interpolater.H>

namespace BoxLib
{
    bool ProperlyNested (const IntVect& ratio, int blockint_factor, int ngrow, 
			 const IndexType& boxType, Interpolater* mapper);

    void FillPatchSingleLevel (MultiFab& mf, Real time, 
			       const PArray<MultiFab>& smf, const std::vector<Real>& stime, 
			       int scomp, int dcomp, int ncomp,
			       const Geometry& geom, PhysBCFunct& physbcf);

    void FillPatchTwoLevels (MultiFab& mf, Real time,
			     const PArray<MultiFab>& cmf, const std::vector<Real>& ct,
			     const PArray<MultiFab>& fmf, const std::vector<Real>& ft,
			     int scomp, int dcomp, int ncomp,
			     const Geometry& cgeom, const Geometry& fgeom, 
			     PhysBCFunct& cbc, PhysBCFunct& fbc,
			     const IntVect& ratio, 
			     Interpolater* mapper, const Array<BCRec>& bcs);

    void InterpFromCoarseLevel (MultiFab& mf, Real time,
				const MultiFab& cmf, int scomp, int dcomp, int ncomp,
				const Geometry& cgeom, const Geometry& fgeom, 
				PhysBCFunct& cbc, PhysBCFunct& fbc, const IntVect& ratio, 
				Interpolater* mapper, const Array<BCRec>& bcs);
}

#endif
