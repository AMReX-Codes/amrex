#ifdef USEHYPRE

#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BndryData.H>
#include <AMReX_MacBndry.H>
#include <AMReX_MultiFabUtil.H>

#include <HypreABecLap.H>

using namespace amrex;

void setBndryConds(BndryData& levelbd, int ibnd, IntVect ratio);

void solve_with_hypre(const Vector<MultiFab*>& soln, Real a, Real b, 
		      const Vector<MultiFab*>& alph, 
		      const Vector<MultiFab*>& beta, 
		      const Vector<MultiFab*>& rhs, 
		      const Vector<Geometry>& geom, 
		      const Vector<BoxArray>& grids,
		      int ibnd)
{
  const Real run_strt = ParallelDescriptor::second();

  int composite_solve = 0;
  Real tolerance_rel, tolerance_abs;
  int max_iter = 100;
  {
    ParmParse pp;
    pp.query("composite_solve", composite_solve);
    
    pp.get("tol_rel", tolerance_rel);
    pp.get("tol_abs", tolerance_abs);

    pp.query("max_iter", max_iter);
  }

  int nlevel = geom.size();

  IntVect ratio = 2.*IntVect::TheUnitVector();

  if (! composite_solve) {
    std::cout << "Non-composite-solve not implemented in Hypre.  Change to composite solve." << std::endl;
    composite_solve = 1;
  }

  if (composite_solve) {
    HypreABecLap hypreSolver(0, nlevel-1, 1);

    // add grids in reverse order in case we are masking covered part of coarse
    for (int level = nlevel-1; level>=0; level--) {
      hypreSolver.addLevel(level, geom[level], grids[level], ratio);
    }

    for (int level = 0; level < nlevel; level++) {
      BndryData *bdp = new BndryData(grids[level], 1, geom[level]);
      setBndryConds(*bdp, ibnd, ratio);

      hypreSolver.setBndry(level, *bdp);
    }

    hypreSolver.buildMatrixStructure();

    hypreSolver.setScalars(a, b);

    for (int level = 0; level < nlevel; level++) {

      hypreSolver.setACoeffs(level, *alph[level]);

      Vector<std::unique_ptr<MultiFab> > bcoeffs(BL_SPACEDIM);
      for (int n = 0; n < BL_SPACEDIM ; n++) 
      {
	  BoxArray edge_boxes(grids[level]);
	  edge_boxes.surroundingNodes(n);
	  bcoeffs[n].reset(new MultiFab(edge_boxes, 1, 0));
      }

      amrex::average_cellcenter_to_face(amrex::GetVecOfPtrs(bcoeffs), *beta[level], geom[level]);
	
      for (int n = 0; n < BL_SPACEDIM ; n++) 
      {
	  hypreSolver.setBCoeffs(level, *bcoeffs[n], n);
      }

      hypreSolver.setRhs(level, *rhs[level]);

      hypreSolver.setInitGuess(level, *soln[level]);
    }

    hypreSolver.solve(soln, tolerance_rel, tolerance_abs, max_iter);
  }
  else {

  }

  Real run_time = ParallelDescriptor::second() - run_strt;

  ParallelDescriptor::ReduceRealMax(run_time, ParallelDescriptor::IOProcessorNumber());
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Total Hypre Run time      : " << run_time << std::endl;
  }
}

void setBndryConds(BndryData& levelbd, int ibnd, IntVect ratio)
{
  int comp = 0;
  Real bc_value = 0.0; // This is hardwired.

  const BoxArray& grids = levelbd.boxes();
  const Geometry& geom = levelbd.getGeom();
  const Real* dx = geom.CellSize();
  const Box& domain = levelbd.getDomain();

  for (OrientationIter fi; fi; ++fi) {
    Orientation face(fi());

    int dir = face.coordDir();
    Real delta = dx[dir]*ratio[dir];
    
    for (FabSetIter bfsi(levelbd[face]); bfsi.isValid(); ++bfsi) {
      const int i = bfsi.index();
      const Box& grd = grids[i];

      if (domain[face] == grd[face] && !geom.isPeriodic(dir)) {
	levelbd.setBoundCond(face, i, comp, ibnd);
	levelbd.setBoundLoc(face, i, 0.0);
	levelbd.setValue(face, i, bc_value);
      }
      else {
	// internal bndry
	levelbd.setBoundCond(face, i, comp, LO_DIRICHLET);
	levelbd.setBoundLoc(face, i, 0.5*delta);
      }
    }
  }
}

#endif
