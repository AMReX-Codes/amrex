#include <Utility.H>
#include <ParmParse.H>
#include <PArray.H>
#include <LO_BCTYPES.H>
#include <MultiFab.H>
#include <Geometry.H>

#include <MacBndry.H>
#include <MGT_Solver.H>
#include <stencil_types.H>

#include <COEF_F.H>

void solve_with_F90(PArray<MultiFab>& soln, Real a, Real b, 
		    const PArray<MultiFab>& alph, 
		    const PArray<MultiFab>& beta, 
		    PArray<MultiFab>& rhs, 
		    const std::vector<Geometry>& geom, 
		    const std::vector<BoxArray>& grids,
		    int ibnd)
{
  const Real run_strt = ParallelDescriptor::second();

  int composite_solve = 0;
  Real tolerance_rel, tolerance_abs;
  {
    ParmParse pp;
    pp.query("composite_solve", composite_solve);
    
    pp.get("tol_rel", tolerance_rel);
    pp.get("tol_abs", tolerance_abs);
  }

  int nlevel = geom.size();
  
  int mg_bc[2*BL_SPACEDIM];
  
  if (ibnd == LO_NEUMANN) {
    for ( int n = 0; n < BL_SPACEDIM; ++n ) {
      mg_bc[2*n + 0] = MGT_BC_NEU;
      mg_bc[2*n + 1] = MGT_BC_NEU;
    }
  }
  else if (ibnd == LO_DIRICHLET) {
    for ( int n = 0; n < BL_SPACEDIM; ++n ) {
      mg_bc[2*n + 0] = MGT_BC_DIR;
      mg_bc[2*n + 1] = MGT_BC_DIR;
    }
  }
  else { // periodic
    for ( int n = 0; n < BL_SPACEDIM; ++n ) {
      mg_bc[2*n + 0] = 0; // MGT_BC_PER;
      mg_bc[2*n + 1] = 0; // MGT_BC_PER;
    }
  }
   
  int phys_bc_type;
  // MacBndry will convert phys_bc_type back to LO_DIRICHLET or LO_NEUMANN
  if (ibnd == LO_DIRICHLET) {
    phys_bc_type = Outflow; // This will be converted to LO_DIRICHLET
  }
  else if (ibnd == LO_NEUMANN) {
    phys_bc_type = Symmetry; // This will be converted to LO_NEUMANN
  }
  else { // periodic
    phys_bc_type = Interior;
  }
  
  BCRec phys_bc;
  for (int n=0; n<BL_SPACEDIM; n++) {
    phys_bc.setLo(n, phys_bc_type);
    phys_bc.setHi(n, phys_bc_type);
  }

  std::vector<DistributionMapping> dmap(nlevel);
  for (int ilev=0; ilev<nlevel; ilev++ ) {
    dmap[ilev] = soln[ilev].DistributionMap();
  }

  // The coefficients are set such that we will solve
  //  (a alpha - b del dot beta grad) soln = rhs
  //  written in the form 
  //  (acoeffs - b del dot bcoeffs grad) soln = rhs

  PArray<MultiFab> acoeffs(nlevel, PArrayManage);
  for (int ilev=0; ilev<nlevel; ilev++ ) {
    acoeffs.set(ilev, new MultiFab(grids[ilev], 1, 0, Fab_allocate));
    acoeffs[ilev].copy(alph[ilev]);
    acoeffs[ilev].mult(a); 
  }

  Array< PArray<MultiFab> > bcoeffs(nlevel);
  for (int ilev=0; ilev<nlevel; ilev++ ) {

    bcoeffs[ilev].resize(BL_SPACEDIM, PArrayManage);

    for (int n = 0; n < BL_SPACEDIM ; n++) {

      BoxArray edge_boxes(grids[ilev]);
      edge_boxes.surroundingNodes(n);
      
      bcoeffs[ilev].set(n, new MultiFab(edge_boxes,1,0,Fab_allocate));

      for (MFIter mfi(bcoeffs[ilev][n]); mfi.isValid(); ++mfi) {
  	int i = mfi.index();
  	const Box& bx = grids[ilev][i];
  	const int* betalo = beta[ilev][mfi].loVect();
  	const int* betahi = beta[ilev][mfi].hiVect();
  	const int* edgelo = bcoeffs[ilev][n][mfi].loVect();
  	const int* edgehi = bcoeffs[ilev][n][mfi].hiVect();
	  
  	FORT_COEF_TO_EDGES(&n, bcoeffs[ilev][n][mfi].dataPtr(),
  			   ARLIM(edgelo), ARLIM(edgehi),
  			   beta[ilev][mfi].dataPtr(),
  			   ARLIM(betalo), ARLIM(betahi),
  			   bx.loVect(),bx.hiVect());
      }
    }
  }

  Array< Array<Real> > xa(nlevel);
  Array< Array<Real> > xb(nlevel);
  for (int ilev=0; ilev<nlevel; ilev++ ) {
    xa[ilev].resize(BL_SPACEDIM);
    xb[ilev].resize(BL_SPACEDIM);
    if (ilev == 0) {
      // For level 0, the boundary lives exactly on the faces
      for (int n=0; n<BL_SPACEDIM; n++) {
	xa[0][n] = 0.0;
	xb[0][n] = 0.0;
      }
    }
    else {
      const Real* dx_crse = geom[ilev-1].CellSize();
      for (int n=0; n<BL_SPACEDIM; n++) {
	xa[ilev][n] = 0.5 * dx_crse[n];
	xb[ilev][n] = 0.5 * dx_crse[n];
      }
    }
  }


  if (composite_solve) {

    bool nodal = false;
    int stencil = CC_CROSS_STENCIL;
    MGT_Solver mgt_solver(geom, mg_bc, grids, dmap, nodal, stencil);

    int index_order = 1; // because of bcoeffs[nlevel][BL_SPACEDIM]
    mgt_solver.set_visc_coefficients(acoeffs, bcoeffs, b, xa, xb, index_order);
    
    MultiFab* soln_p[nlevel];
    MultiFab* rhs_p[nlevel];
    for (int ilev=0; ilev<nlevel; ilev++ ) {
      soln_p[ilev] = &soln[ilev];
      rhs_p[ilev] = &rhs[ilev];
    }    

    MacBndry bndry(grids[0], 1, geom[0]);
    bndry.setBndryValues(soln[0], 0, 0, 1, phys_bc); 
    // does this work for Neumann?

    int always_use_bnorm = 0;
    Real final_resnorm;
    mgt_solver.solve(soln_p, rhs_p, bndry, tolerance_rel, tolerance_abs, always_use_bnorm, final_resnorm);
 
  }
  else {

    for (int ilev=0; ilev<nlevel; ilev++ ) {

      std::vector<Geometry> geom_l(1, geom[ilev]);
      std::vector<BoxArray> grids_l(1, grids[ilev]);
      std::vector<DistributionMapping> dmap_l(1, dmap[ilev]);

      bool nodal = false;
      int stencil = CC_CROSS_STENCIL;
      MGT_Solver mgt_solver(geom_l, mg_bc, grids_l, dmap_l, nodal, stencil);

      MultiFab* acoeffs_l[1];
      acoeffs_l[0] = &acoeffs[ilev];

      MultiFab* bcoeffs_l[1][BL_SPACEDIM];
      for (int n = 0; n < BL_SPACEDIM ; n++) {
	bcoeffs_l[0][n] = &bcoeffs[ilev][n];
      }
      
      Array< Array<Real> > xa_l(1, xa[ilev]);
      Array< Array<Real> > xb_l(1, xb[ilev]);

      mgt_solver.set_visc_coefficients(acoeffs_l, bcoeffs_l, b, xa_l, xb_l);

      MultiFab* soln_p[1];
      MultiFab* rhs_p[1];
      soln_p[0] = &soln[ilev];
      rhs_p[0] = &rhs[ilev];

      MacBndry bndry(grids[ilev], 1, geom[ilev]);

      const int src_comp  = 0;
      const int dest_comp = 0;
      const int num_comp  = 1;

      if (ilev == 0) {
	bndry.setBndryValues(soln[ilev], src_comp, dest_comp, num_comp, phys_bc); 
	// does this work for Neumann?
      }
      else {
	IntVect crse_ratio = 2.*IntVect::TheUnitVector();
	const int num_comp = 1;
	const int in_rad     = 0;
        const int out_rad    = 1;
        const int extent_rad = 2;
	BoxArray crse_boxes = BoxArray(grids[ilev]).coarsen(crse_ratio);

	BndryRegister crse_br(crse_boxes, in_rad, out_rad, extent_rad, num_comp);
        crse_br.copyFrom(soln[ilev-1], extent_rad, src_comp, dest_comp, num_comp);

	bndry.setBndryValues(crse_br, src_comp, soln[ilev], src_comp,
			     dest_comp, num_comp, crse_ratio, phys_bc);
      }

      Real final_resnorm;
      int always_use_bnorm = 0;
      mgt_solver.solve(soln_p, rhs_p, bndry, tolerance_rel, tolerance_abs, always_use_bnorm, final_resnorm);
    }
  }

  Real run_time = ParallelDescriptor::second() - run_strt;

  ParallelDescriptor::ReduceRealMax(run_time, ParallelDescriptor::IOProcessorNumber());
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Total BoxLib_F Run time      : " << run_time << std::endl;
  }
}

