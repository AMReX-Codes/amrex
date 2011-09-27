#include <Utility.H>
#include <ParmParse.H>
#include <PArray.H>
#include <LO_BCTYPES.H>
#include <MultiFab.H>
#include <Geometry.H>

#include <MacBndry.H>
#include <MGT_Solver.H>

#include "COEF_F.H"

void solve_with_F90(PArray<MultiFab>& soln, int iF90, Real a, Real b, 
		    const PArray<MultiFab>& alph, 
		    const PArray<MultiFab>& beta, 
		    PArray<MultiFab>& rhs, 
		    const std::vector<Geometry>& geom, 
		    const std::vector<BoxArray>& grids,
		    int ibnd)
{
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
      mg_bc[2*n + 0] = MGT_BC_PER;
      mg_bc[2*n + 1] = MGT_BC_PER;
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
 
  if (composite_solve) {

  }
  else {

    for (int ilev=0; ilev<nlevel; ilev++ ) {

      std::vector<Geometry> fgeom(1);
      fgeom[0] = geom[ilev];

      std::vector<BoxArray> fgrids(1);
      fgrids[0] = grids[ilev];

      std::vector<DistributionMapping> dmap(1);
      dmap[0] = soln[ilev].DistributionMap();

      bool nodal = false;
      MGT_Solver mgt_solver(fgeom, mg_bc, fgrids, dmap, nodal);

      PArray<MultiFab> acoeffs(1, PArrayManage);
      acoeffs.set(0, new MultiFab(grids[ilev], 1, 0, Fab_allocate));
      acoeffs[0].copy(alph[ilev]);
      acoeffs[0].mult(a); 

      Array< PArray<MultiFab> > bcoeffs(BL_SPACEDIM);
      for (int n = 0; n < BL_SPACEDIM ; n++) {

	bcoeffs[n].resize(1, PArrayManage);

	BoxArray edge_boxes(grids[ilev]);
	edge_boxes.surroundingNodes(n);
      
	bcoeffs[n].set(0, new MultiFab(edge_boxes,1,0,Fab_allocate));

	for (MFIter mfi(bcoeffs[n][0]); mfi.isValid(); ++mfi) {
	  int i = mfi.index();
	  const Box& bx = grids[ilev][i];
	  const int* betalo = beta[ilev][i].loVect();
	  const int* betahi = beta[ilev][i].hiVect();
	  const int* edgelo = bcoeffs[n][0][i].loVect();
	  const int* edgehi = bcoeffs[n][0][i].hiVect();
	  
	  FORT_COEF_TO_EDGES(&n, bcoeffs[n][0][i].dataPtr(),
			     ARLIM(edgelo), ARLIM(edgehi),
			     beta[ilev][i].dataPtr(),
			     ARLIM(betalo), ARLIM(betahi),
			     bx.loVect(),bx.hiVect());
	}
      }

      Array< Array<Real> > xa(1);
      Array< Array<Real> > xb(1);
      xa[0].resize(BL_SPACEDIM);
      xb[0].resize(BL_SPACEDIM);
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
	  xa[0][n] = 0.5 * dx_crse[n];
	  xb[0][n] = 0.5 * dx_crse[n];
	}
      }

      // The coefficients are set such that we will solve
      //  (a alpha - b del dot beta grad) soln = rhs
      //  written in the form 
      //  (acoeffs - b del dot bcoeffs grad) soln = rhs
      mgt_solver.set_visc_coefficients(acoeffs, bcoeffs, b, xa, xb);

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
      mgt_solver.solve(soln_p, rhs_p, tolerance_rel, tolerance_abs, bndry, final_resnorm);
    }
  }
}

