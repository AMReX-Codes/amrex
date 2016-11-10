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
			       const Array<MultiFab*>& smf, const Array<Real>& stime, 
			       int scomp, int dcomp, int ncomp,
			       const Geometry& geom, PhysBCFunctBase& physbcf);

    void FillPatchTwoLevels (MultiFab& mf, Real time,
			     const Array<MultiFab*>& cmf, const Array<Real>& ct,
			     const Array<MultiFab*>& fmf, const Array<Real>& ft,
			     int scomp, int dcomp, int ncomp,
			     const Geometry& cgeom, const Geometry& fgeom, 
			     PhysBCFunctBase& cbc, PhysBCFunctBase& fbc,
			     const IntVect& ratio, 
			     Interpolater* mapper, const Array<BCRec>& bcs);

    void InterpFromCoarseLevel (MultiFab& mf, Real time,
				const MultiFab& cmf, int scomp, int dcomp, int ncomp,
				const Geometry& cgeom, const Geometry& fgeom, 
				PhysBCFunctBase& cbc, PhysBCFunctBase& fbc, const IntVect& ratio, 
				Interpolater* mapper, const Array<BCRec>& bcs);
}

#endif
