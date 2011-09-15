// We solve
// (a alpha - b del dot beta grad) soln = rhs
// where a and b are scalars, alpha and beta are arrays

#include <fstream>
#include <iomanip>

#include <Utility.H>
#include <ParmParse.H>
#include <LO_BCTYPES.H>
#include <BndryData.H>
#include <MultiGrid.H>
#include <CGSolver.H>
#include <Laplacian.H>
#include <ABecLaplacian.H>
#include <ParallelDescriptor.H>
#include <VisMF.H>
#include <COEF_F.H>
#include <MacBndry.H>
#include <MGT_Solver.H>
#ifdef USEHYPRE
#include <HypreABecLap.H>
#endif

int  verbose       = 2;     
Real tolerance_rel = 1.e-12;
Real tolerance_abs = 1.e-12;
int  maxiter       = 100; 
int  plot_rhs      = 0; 
int  plot_soln     = 0; 
int  plot_asol     = 0; 
int  plot_err      = 0;
int  comp_norm     = 1;

Real dx[BL_SPACEDIM];

std::string solver_type = "fillme";
std::string     bc_type = "fillme";

int Ncomp = 1;

void compute_analyticSolution(MultiFab& anaSoln)
{
  int ibnd;
  if (bc_type == "Dirichlet") {
    ibnd = LO_DIRICHLET;
  }
  else if (bc_type == "Neumann") {
    ibnd = LO_NEUMANN;
  }
  else if (bc_type == "Periodic") {
    ibnd = 0;
  }

  for (MFIter mfi(anaSoln); mfi.isValid(); ++mfi) {
    const int* alo = anaSoln[mfi].loVect();
    const int* ahi = anaSoln[mfi].hiVect();
    const Box& bx = mfi.validbox();

    FORT_COMP_ASOL(anaSoln[mfi].dataPtr(), ARLIM(alo), ARLIM(ahi),
		   bx.loVect(),bx.hiVect(),dx, &ibnd);
  }
}

void setup_coeffs(BoxArray& bs, MultiFab& alpha, MultiFab beta[])
{
  MultiFab cc_coef(bs,Ncomp,1); // cell-centered beta

  for ( MFIter mfi(alpha); mfi.isValid(); ++mfi ) {
    const Box& bx = mfi.validbox();

    const int* alo = alpha[mfi].loVect();
    const int* ahi = alpha[mfi].hiVect();
    FORT_SET_ALPHA(alpha[mfi].dataPtr(),ARLIM(alo),ARLIM(ahi),
    		   bx.loVect(),bx.hiVect(),dx);

    const int* clo = cc_coef[mfi].loVect();
    const int* chi = cc_coef[mfi].hiVect();
    FORT_SET_CC_COEF(cc_coef[mfi].dataPtr(),ARLIM(clo),ARLIM(chi),
		     bx.loVect(),bx.hiVect(),dx);
  }

  // convert cell-centered beta to edges
  for ( int n=0; n<BL_SPACEDIM; ++n ) {
    for ( MFIter mfi(beta[n]); mfi.isValid(); ++mfi ) {
      int i = mfi.index();
      Box bx(bs[i]);
      const int* clo = cc_coef[i].loVect();
      const int* chi = cc_coef[i].hiVect();
      const int* edgelo = beta[n][i].loVect();
      const int* edgehi = beta[n][i].hiVect();
      
      FORT_COEF_TO_EDGES(&n,beta[n][mfi].dataPtr(),ARLIM(edgelo),ARLIM(edgehi),
			 cc_coef[mfi].dataPtr(),ARLIM(clo),ARLIM(chi),
			 bx.loVect(),bx.hiVect());
    }
  }
}

void setup_rhs(MultiFab& rhs)
{
  int ibnd;
  if (bc_type == "Dirichlet") {
    ibnd = LO_DIRICHLET;
  }
  else if (bc_type == "Neumann") {
    ibnd = LO_NEUMANN;
  }
  else if (bc_type == "Periodic") {
    ibnd = 0;
  }

  for ( MFIter mfi(rhs); mfi.isValid(); ++mfi ) {
    const int* rlo = rhs[mfi].loVect();
    const int* rhi = rhs[mfi].hiVect();
    const Box& bx = mfi.validbox();

    FORT_SET_RHS(rhs[mfi].dataPtr(),ARLIM(rlo),ARLIM(rhi),
                 bx.loVect(),bx.hiVect(),dx, &ibnd);
  }

  if (plot_rhs == 1) {
    VisMF::Write(rhs,"RHS");
  }
}

void solve_with_Cpp(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab beta[], 
		    MultiFab& rhs, const BoxArray& bs, const Geometry& geom)
{
  BndryData bd(bs, 1, geom);
  int comp = 0;

  for (int n=0; n<BL_SPACEDIM; ++n) {
    for ( MFIter mfi(rhs); mfi.isValid(); ++mfi ) {
      int i = mfi.index(); 

      // Set the boundary conditions to live exactly on the faces of the domain
      bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.0 );
      bd.setBoundLoc(Orientation(n, Orientation::high),i,0.0 );

      if (bc_type == "Dirichlet") {
	// Define the type of boundary conditions to be Dirichlet
	bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,LO_DIRICHLET);
	bd.setBoundCond(Orientation(n, Orientation::high),i,comp,LO_DIRICHLET);
      }
      else if (bc_type == "Neumann") {
	// Define the type of boundary conditions to be Neumann
	bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,LO_NEUMANN);
	bd.setBoundCond(Orientation(n, Orientation::high),i,comp,LO_NEUMANN);
      }

      // Set the Dirichlet/Neumann values to be 0
      bd.setValue(Orientation(n, Orientation::low) ,i, 0.0);
      bd.setValue(Orientation(n, Orientation::high),i, 0.0);
    }
  }

  ABecLaplacian abec_operator(bd, dx);
  abec_operator.setScalars(a, b);
  abec_operator.setCoefficients(alpha, beta);

  MultiGrid mg(abec_operator);
  mg.setVerbose(verbose);
  mg.solve(soln, rhs, tolerance_rel, tolerance_abs);
} 

void solve_with_F90(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab beta[], 
		    MultiFab& rhs, const BoxArray& bs, const Geometry& geom)
{
  // Translate into F90 solver
  std::vector<BoxArray> bav(1);
  bav[0] = bs;
  std::vector<DistributionMapping> dmv(1);
  dmv[0] = rhs.DistributionMap();
  std::vector<Geometry> fgeom(1);
  fgeom[0] = geom;

  int mg_bc[2*BL_SPACEDIM];

  Array< Array<Real> > xa(1);
  Array< Array<Real> > xb(1);

  xa[0].resize(BL_SPACEDIM);
  xb[0].resize(BL_SPACEDIM);

  // Set the boundary conditions to live exactly on the faces of the domain
  for ( int dir = 0; dir < BL_SPACEDIM; ++dir ) {
    xa[0][dir] = 0.;
    xb[0][dir] = 0.;
  }

  if (bc_type == "Periodic") {
    // Define the type of boundary conditions to be periodic
    for ( int dir = 0; dir < BL_SPACEDIM; ++dir ) {
      mg_bc[2*dir + 0] = MGT_BC_PER;
      mg_bc[2*dir + 1] = MGT_BC_PER;
    }
  }
  else if (bc_type == "Neumann") {
    // Define the type of boundary conditions to be Neumann
    for ( int dir = 0; dir < BL_SPACEDIM; ++dir ) {
      mg_bc[2*dir + 0] = MGT_BC_NEU;
      mg_bc[2*dir + 1] = MGT_BC_NEU;
    }
  }
  else if (bc_type == "Dirichlet") {
    // Define the type of boundary conditions to be Dirichlet
    for ( int dir = 0; dir < BL_SPACEDIM; ++dir ) {
      mg_bc[2*dir + 0] = MGT_BC_DIR;
      mg_bc[2*dir + 1] = MGT_BC_DIR;
    }
  }

  MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false);

  MultiFab* soln_p[1]; soln_p[0] = &soln;
  MultiFab*  rhs_p[1];  rhs_p[0] = &rhs;

  PArray<MultiFab> acoeffs(1, PArrayManage);
  acoeffs.set(0, new MultiFab(bs,Ncomp,0,Fab_allocate));
  acoeffs[0].copy(alpha);
  acoeffs[0].mult(a); 

  Array< PArray<MultiFab> > bcoeffs(BL_SPACEDIM);
  for (int i = 0; i < BL_SPACEDIM ; i++) {
    BoxArray edge_boxes(bs);
    edge_boxes.surroundingNodes(i);

    bcoeffs[i].resize(1,PArrayManage);
    bcoeffs[i].set(0, new MultiFab(edge_boxes,Ncomp,0,Fab_allocate));
    bcoeffs[i][0].copy(beta[i]);
  }

  // The coefficients are set such that we will solve
  //  (a alpha - b del dot beta grad) soln = rhs
  //  written in the form 
  //  (acoeffs - b del dot bcoeffs grad) soln = rhs
  mgt_solver.set_visc_coefficients(acoeffs,bcoeffs,b,xa,xb);

  BCRec* phys_bc = new BCRec;
  for (int i = 0; i < BL_SPACEDIM; i++) {
    phys_bc->setLo(i,0);
    phys_bc->setHi(i,0);
  }

  MacBndry bndry(bs,1,geom);
  bndry.setBndryValues(soln,0,0,1,*phys_bc);
  
  Real final_resnorm;
  mgt_solver.solve(soln_p, rhs_p, tolerance_rel, tolerance_abs, bndry, final_resnorm);
}

#ifdef USEHYPRE
void solve_with_hypre(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab beta[], 
		      MultiFab& rhs, const BoxArray& bs, const Geometry& geom)
{
  BndryData bd(bs, 1, geom);
  int comp = 0;

  for (int n=0; n<BL_SPACEDIM; ++n) {
    for ( MFIter mfi(rhs); mfi.isValid(); ++mfi ) {
      int i = mfi.index(); 
      
      // Set the boundary conditions to live exactly on the faces of the domain
      bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.0 );
      bd.setBoundLoc(Orientation(n, Orientation::high),i,0.0 );
      
      if (bc_type == "Dirichlet") {
	// Define the type of boundary conditions to be Dirichlet
	bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,LO_DIRICHLET);
	bd.setBoundCond(Orientation(n, Orientation::high),i,comp,LO_DIRICHLET);
      }
      else if (bc_type == "Neumann") {
	// Define the type of boundary conditions to be Neumann
	bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,LO_NEUMANN);
	bd.setBoundCond(Orientation(n, Orientation::high),i,comp,LO_NEUMANN);
      }
      
      // Set the Dirichlet/Neumann values to be 0
      bd.setValue(Orientation(n, Orientation::low) ,i, 0.0);
      bd.setValue(Orientation(n, Orientation::high),i, 0.0);
    }
  }

  HypreABecLap hypreSolver(bs, geom);
  hypreSolver.setScalars(a, b);
  hypreSolver.setACoeffs(alpha);
  hypreSolver.setBCoeffs(beta);
  hypreSolver.setVerbose(verbose);
  hypreSolver.solve(soln, rhs, tolerance_rel, tolerance_abs, maxiter, bd);
}
#endif

int main(int argc, char* argv[])
{
  BoxLib::Initialize(argc,argv);

  std::cout << std::setprecision(15);

  ParmParse ppmg("mg");  
  ppmg.query("v", verbose);
  
  ParmParse pp;
  pp.query("tol_rel", tolerance_rel);
  pp.query("tol_abs", tolerance_abs);
  pp.query("maxiter", maxiter);
  pp.query("plot_rhs" , plot_rhs);
  pp.query("plot_soln", plot_soln);
  pp.query("plot_asol", plot_asol);
  pp.query("plot_err", plot_err);
  pp.query("comp_norm", comp_norm);

  Real a = 0., b = 1.;
  pp.query("a",  a);
  pp.query("b",  b);

  int n_cell;
  int max_grid_size;

  pp.get("n_cell",n_cell);
  pp.get("max_grid_size",max_grid_size);
  pp.get("solver_type",solver_type);
  pp.get("bc_type",bc_type);

  // Define a single box covering the domain
  IntVect dom_lo(0,0,0);
  IntVect dom_hi(n_cell-1,n_cell-1,n_cell-1);
  Box domain(dom_lo,dom_hi);

  // Initialize the boxarray "bs" from the single box "bx"
  BoxArray bs(domain);

  // Break up boxarray "bs" into chunks no larger than "max_grid_size" along a direction
  bs.maxSize(max_grid_size);

  // This defines the physical size of the box.  Right now the box is [0,1] in each direction.
  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    real_box.setLo(n, 0.0);
    real_box.setHi(n, 1.0);
  }
 
  // This says we are using Cartesian coordinates
  int coord = 0;
  
  // This sets the boundary conditions to be periodic or not
  int* is_per = new int[BL_SPACEDIM];
  
  if (bc_type == "Dirichlet" || bc_type == "Neumann") {
    for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 0;
  } 
  else if (bc_type == "Periodic") {
    for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 1;
  }
  else {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Don't know this BC type: " << bc_type << std::endl;
    }
    BoxLib::Error("");
  }
 
  // This defines a Geometry object which is useful for writing the plotfiles
  Geometry geom(domain,&real_box,coord,is_per);

  for ( int n=0; n<BL_SPACEDIM; n++ ) {
    dx[n] = ( geom.ProbHi(n) - geom.ProbLo(n) )/domain.length(n);
  }

  // Allocate and define the right hand side.
  MultiFab rhs(bs, Ncomp, 0, Fab_allocate); 
  setup_rhs(rhs);

  // Allocate the solution array and initialize to zero
  // Set the number of ghost cells in the solution array.
  MultiFab soln(bs, Ncomp, 1, Fab_allocate);
  soln.setVal(0.0);

  MultiFab alpha(bs, Ncomp, 0, Fab_allocate);
  MultiFab beta[BL_SPACEDIM];
  for ( int n=0; n<BL_SPACEDIM; ++n ) {
    BoxArray bx(bs);
    beta[n].define(bx.surroundingNodes(n), Ncomp, 1, Fab_allocate);
  }

  setup_coeffs(bs, alpha, beta);

  const Real run_strt = ParallelDescriptor::second();

  if (solver_type == "BoxLib_C") {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Solving with BoxLib C++ solver " << std::endl;
    }
    solve_with_Cpp(soln, a, b, alpha, beta, rhs, bs, geom);
  } 
  else if (solver_type == "BoxLib_F") {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Solving with BoxLib Fortran90 solver " << std::endl;
    }
    solve_with_F90(soln, a, b, alpha, beta, rhs, bs, geom);
  }
  else if (solver_type == "Hypre") {
#ifdef USEHYPRE
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Solving with Hypre " << std::endl;
    }
    solve_with_hypre(soln, a, b, alpha, beta, rhs, bs, geom);
#else
    BoxLib::Error("Set USE_HYPRE=TRUE in GNUMakefile.");
#endif
  }
  else {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Don't know this solver type: " << solver_type << std::endl;
    }
    BoxLib::Error("");
  }

  Real run_stop = ParallelDescriptor::second() - run_strt;
  
  ParallelDescriptor::ReduceRealMax(run_stop,ParallelDescriptor::IOProcessorNumber());
  
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Run time = " << run_stop << std::endl;
  }
  
  if (plot_soln == 1) {
    VisMF::Write(soln,"SOLN");
  }

  if (comp_norm || plot_err || plot_asol) {
    MultiFab anaSoln(bs, Ncomp, 0, Fab_allocate);
    compute_analyticSolution(anaSoln);
    
    if (plot_asol) {
      VisMF::Write(anaSoln,"ASOL");
    }

    anaSoln.minus(soln, 0, Ncomp, 0); // anaSoln is error now
    if (plot_err) {
      VisMF::Write(anaSoln,"ERR");
    }

    if (comp_norm) {
      Real twoNorm = anaSoln.norm2();
      Real maxNorm = anaSoln.norm0();

      anaSoln.setVal(1.0);
      Real vol = anaSoln.norm2();
      twoNorm /= vol;

      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "2 norm error: " << twoNorm << ", max norm error: " << maxNorm << std::endl;
      }
    }
  }

  BoxLib::Finalize();
}
