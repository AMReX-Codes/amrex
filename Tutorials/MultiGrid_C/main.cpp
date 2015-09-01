// We solve (a alpha - b del dot beta grad) soln = rhs
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
#include <MacBndry.H>
#ifdef USE_F90_SOLVERS
#include <MGT_Solver.H>
#include <stencil_types.H>
#endif

#ifdef USEHYPRE
#include <HypreABecLap.H>
#endif

#include <COEF_F.H>
#include <RHS_F.H>
#include <writePlotFile.H>

int  verbose       = 2;     
Real tolerance_rel = 1.e-8;
Real tolerance_abs = 0.0;
int  maxiter       = 100; 
int  plot_rhs      = 0; 
int  plot_beta     = 0;
int  plot_soln     = 0; 
int  plot_asol     = 0; 
int  plot_err      = 0;
int  comp_norm     = 1;

Real dx[BL_SPACEDIM];

enum solver_t {BoxLib_C, BoxLib_F, Hypre, HPGMG, All};
enum bc_t {Periodic = 0,
	   Dirichlet = LO_DIRICHLET, 
	   Neumann = LO_NEUMANN};
solver_t solver_type;
bc_t     bc_type;

int Ncomp = 1;
int n_cell;
int max_grid_size;

int domain_boundary_condition;

void compute_analyticSolution(MultiFab& anaSoln);
void setup_coeffs(BoxArray& bs, MultiFab& alpha, MultiFab beta[], const Geometry& geom, MultiFab& beta_cc);
void setup_rhs(MultiFab& rhs, const Geometry& geom);
void set_boundary(BndryData& bd, const MultiFab& rhs, int comp);
void solve(MultiFab& soln, const MultiFab& anaSoln, 
	   Real a, Real b, MultiFab& alpha, MultiFab beta[], MultiFab& beta_cc,
	   MultiFab& rhs, const BoxArray& bs, const Geometry& geom,
	   solver_t solver);
void solve_with_Cpp(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab beta[], 
		    MultiFab& rhs, const BoxArray& bs, const Geometry& geom);

#ifdef USE_F90_SOLVERS
void solve_with_F90(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab beta[], 
		    MultiFab& rhs, const BoxArray& bs, const Geometry& geom);
#endif

#ifdef USEHYPRE
void solve_with_hypre(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab beta[], 
		      MultiFab& rhs, const BoxArray& bs, const Geometry& geom);
#endif

#ifdef USEHPGMG
void solve_with_HPGMG(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab& beta_cc,
		    MultiFab& rhs, const BoxArray& bs, const Geometry& geom, int n_cell);
#endif

int main(int argc, char* argv[])
{
  BoxLib::Initialize(argc,argv);

  BL_PROFILE_VAR("main()", pmain);

  std::cout << std::setprecision(15);

  ParmParse ppmg("mg");  
  ppmg.query("v", verbose);
  
  ParmParse pp;

  {
    std::string solver_type_s;
    pp.get("solver_type",solver_type_s);
    if (solver_type_s == "BoxLib_C") {
      solver_type = BoxLib_C;
    }
    else if (solver_type_s == "BoxLib_F") {
#ifdef USE_F90_SOLVERS
      solver_type = BoxLib_F;      
#else
      BoxLib::Error("Set USE_FORTRAN=TRUE in GNUmakefile");
#endif
    }
    else if (solver_type_s == "Hypre") {
#ifdef USEHYPRE
      solver_type = Hypre;
#else
      BoxLib::Error("Set USE_HYPRE=TRUE in GNUmakefile");
#endif
    }
    else if (solver_type_s == "All") {
      solver_type = All;
    }  
    else {
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Don't know this solver type: " << solver_type << std::endl;
      }
      BoxLib::Error("");
    }
  }

  {
    std::string bc_type_s;
    pp.get("bc_type",bc_type_s);
    if (bc_type_s == "Dirichlet") {
      bc_type = Dirichlet;
#ifdef USEHPGMG
      domain_boundary_condition = BC_DIRICHLET;
#endif
    }
    else if (bc_type_s == "Neumann") {
      bc_type = Neumann;
#ifdef USEHPGMG
      BoxLib::Error("HPGMG does not support Neumann boundary conditions");
#endif
    }
    else if (bc_type_s == "Periodic") {
      bc_type = Periodic;
#ifdef USEHPGMG
      domain_boundary_condition = BC_PERIODIC;
#endif
    }
    else {
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Don't know this boundary type: " << bc_type << std::endl;
      }
      BoxLib::Error("");
    }
  }

  pp.query("tol_rel", tolerance_rel);
  pp.query("tol_abs", tolerance_abs);
  pp.query("maxiter", maxiter);
  pp.query("plot_rhs" , plot_rhs);
  pp.query("plot_beta", plot_beta);
  pp.query("plot_soln", plot_soln);
  pp.query("plot_asol", plot_asol);
  pp.query("plot_err", plot_err);
  pp.query("comp_norm", comp_norm);

  Real a, b;
  pp.get("a",  a);
  pp.get("b",  b);

  pp.get("n_cell",n_cell);
  pp.get("max_grid_size",max_grid_size);

  // Define a single box covering the domain
  IntVect dom_lo(D_DECL(0,0,0));
  IntVect dom_hi(D_DECL(n_cell-1,n_cell-1,n_cell-1));
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
  int is_per[BL_SPACEDIM];
  
  if (bc_type == Dirichlet || bc_type == Neumann) {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Using Dirichlet or Neumann boundary conditions." << std::endl;
    }
    for (int n = 0; n < BL_SPACEDIM; n++) is_per[n] = 0;
  } 
  else {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Using periodic boundary conditions." << std::endl;
    }
    for (int n = 0; n < BL_SPACEDIM; n++) is_per[n] = 1;
  }
 
  // This defines a Geometry object which is useful for writing the plotfiles
  Geometry geom(domain,&real_box,coord,is_per);

  for ( int n=0; n<BL_SPACEDIM; n++ ) {
    dx[n] = ( geom.ProbHi(n) - geom.ProbLo(n) )/domain.length(n);
  }

  if (ParallelDescriptor::IOProcessor()) {
     std::cout << "Domain size     : " << n_cell << std::endl;
     std::cout << "Max_grid_size   : " << max_grid_size << std::endl;
     std::cout << "Number of grids : " << bs.size() << std::endl;
  }

  // Allocate and define the right hand side.
  MultiFab rhs(bs, Ncomp, 0, Fab_allocate); 
  setup_rhs(rhs, geom);

  // Set up the Helmholtz operator coefficients.
  MultiFab alpha(bs, Ncomp, 0, Fab_allocate);
  MultiFab beta[BL_SPACEDIM];
  for ( int n=0; n<BL_SPACEDIM; ++n ) {
    BoxArray bx(bs);
    beta[n].define(bx.surroundingNodes(n), Ncomp, 1, Fab_allocate);
  }

  // The way HPGMG stores face-centered data is completely different than the
  // way BoxLib does it, and translating between the two directly via indexing
  // magic is a nightmare. Happily, the way this tutorial calculates
  // face-centered values is by first calculating cell-centered values and then
  // interpolating to the cell faces. HPGMG can do the same thing, so rather
  // than converting directly from BoxLib's face-centered data to HPGMG's, just
  // give HPGMG the cell-centered data and let it interpolate itself.

  MultiFab beta_cc(bs,Ncomp,1); // cell-centered beta
  setup_coeffs(bs, alpha, beta, geom, beta_cc);

  MultiFab anaSoln;
  if (comp_norm || plot_err || plot_asol) {
    anaSoln.define(bs, Ncomp, 0, Fab_allocate);
    compute_analyticSolution(anaSoln);
    
    if (plot_asol) {
      writePlotFile("ASOL", anaSoln, geom);
    }
  }

  // Allocate the solution array 
  // Set the number of ghost cells in the solution array.
  MultiFab soln(bs, Ncomp, 1, Fab_allocate);

#ifdef USEHYPRE
  if (solver_type == Hypre || solver_type == All) {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Solving with Hypre " << std::endl;
    }

    solve(soln, anaSoln, a, b, alpha, beta, beta_cc, rhs, bs, geom, Hypre);
  }
#endif

  if (solver_type == BoxLib_C || solver_type == All) {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Solving with BoxLib C++ solver " << std::endl;
    }

    solve(soln, anaSoln, a, b, alpha, beta, beta_cc, rhs, bs, geom, BoxLib_C);
  }

#ifdef USE_F90_SOLVERS
  if (solver_type == BoxLib_F || solver_type == All) {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Solving with BoxLib F90 solver " << std::endl;
    }

    solve(soln, anaSoln, a, b, alpha, beta, beta_cc, rhs, bs, geom, BoxLib_F);
  }
#endif

#ifdef USEHPGMG
  if (solver_type == HPGMG || solver_type == All) {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Solving with HPGMG solver " << std::endl;
    }

    solve(soln, anaSoln, a, b, alpha, beta, beta_cc, rhs, bs, geom, HPGMG);
  }
#endif

  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "----------------------------------------" << std::endl;
  }
  
  BL_PROFILE_VAR_STOP(pmain);

  BoxLib::Finalize();
}

void compute_analyticSolution(MultiFab& anaSoln)
{
  BL_PROFILE("compute_analyticSolution()");
  int ibnd = static_cast<int>(bc_type); 

  for (MFIter mfi(anaSoln); mfi.isValid(); ++mfi) {
    const int* alo = anaSoln[mfi].loVect();
    const int* ahi = anaSoln[mfi].hiVect();
    const Box& bx = mfi.validbox();

    FORT_COMP_ASOL(anaSoln[mfi].dataPtr(), ARLIM(alo), ARLIM(ahi),
		   bx.loVect(),bx.hiVect(),dx, ibnd);
  }
}

void setup_coeffs(BoxArray& bs, MultiFab& alpha, MultiFab beta[], const Geometry& geom, MultiFab& cc_coef)
{
  BL_PROFILE("setup_coeffs()");
  ParmParse pp;

  Real sigma, w;
  pp.get("sigma", sigma);
  pp.get("w", w);

  for ( MFIter mfi(alpha); mfi.isValid(); ++mfi ) {
    const Box& bx = mfi.validbox();

    const int* alo = alpha[mfi].loVect();
    const int* ahi = alpha[mfi].hiVect();
    FORT_SET_ALPHA(alpha[mfi].dataPtr(),ARLIM(alo),ARLIM(ahi),
    		   bx.loVect(),bx.hiVect(),dx);

    const int* clo = cc_coef[mfi].loVect();
    const int* chi = cc_coef[mfi].hiVect();
    FORT_SET_CC_COEF(cc_coef[mfi].dataPtr(),ARLIM(clo),ARLIM(chi),
		     bx.loVect(),bx.hiVect(),dx, sigma, w);
  }

  // convert cell-centered beta to edges
  for ( int n=0; n<BL_SPACEDIM; ++n ) {
    for ( MFIter mfi(beta[n]); mfi.isValid(); ++mfi ) {
      int i = mfi.index();
      Box bx(bs[i]);
      const int* clo = cc_coef[mfi].loVect();
      const int* chi = cc_coef[mfi].hiVect();
      const int* edgelo = beta[n][mfi].loVect();
      const int* edgehi = beta[n][mfi].hiVect();
      
      FORT_COEF_TO_EDGES(&n,beta[n][mfi].dataPtr(),ARLIM(edgelo),ARLIM(edgehi),
			 cc_coef[mfi].dataPtr(),ARLIM(clo),ARLIM(chi),
			 bx.loVect(),bx.hiVect());
    }
  }

  if (plot_beta == 1) {
    writePlotFile("BETA", cc_coef, geom);
  }
}

void setup_rhs(MultiFab& rhs, const Geometry& geom)
{
  BL_PROFILE("setup_rhs()");
  ParmParse pp;

  Real a, b, sigma, w;
  pp.get("a",  a);
  pp.get("b",  b);
  pp.get("sigma", sigma);
  pp.get("w", w);

  int ibnd = static_cast<int>(bc_type); 

  // We test the sum of the RHS to check solvability
  Real sum_rhs = 0.;

  for ( MFIter mfi(rhs); mfi.isValid(); ++mfi ) {
    const int* rlo = rhs[mfi].loVect();
    const int* rhi = rhs[mfi].hiVect();
    const Box& bx = mfi.validbox();

    FORT_SET_RHS(rhs[mfi].dataPtr(),ARLIM(rlo),ARLIM(rhi),
                 bx.loVect(),bx.hiVect(),dx, a, b, sigma, w, ibnd);
    sum_rhs += rhs[mfi].sum(0,1);
  }

  for (int n=0; n < BL_SPACEDIM; n++) {
    sum_rhs *= dx[n]; 
  }
  
  ParallelDescriptor::ReduceRealSum(sum_rhs, ParallelDescriptor::IOProcessorNumber());

  if (ParallelDescriptor::IOProcessor()) {
     std::cout << "Sum of RHS      : " << sum_rhs << std::endl;
  }

  if (plot_rhs == 1) {
    writePlotFile("RHS", rhs, geom);
  }
}

void set_boundary(BndryData& bd, const MultiFab& rhs, int comp=0)
{
  BL_PROFILE("set_boundary()");
  Real bc_value = 0.0;

  for (int n=0; n<BL_SPACEDIM; ++n) {
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi ) {
      int i = mfi.index(); 
      
      const Box& bx = mfi.validbox();
      
      // Our default will be that the face of this grid is either touching another grid
      //  across an interior boundary or a periodic boundary.  We will test for the other
      //  cases below.
      {
	// Define the type of boundary conditions to be Dirichlet (even for periodic)
	bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,LO_DIRICHLET);
	bd.setBoundCond(Orientation(n, Orientation::high),i,comp,LO_DIRICHLET);
	
	// Set the boundary conditions to the cell centers outside the domain
	bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.5*dx[n]);
	bd.setBoundLoc(Orientation(n, Orientation::high),i,0.5*dx[n]);
      }

      // Now test to see if we should override the above with Dirichlet or Neumann physical bc's
      if (bc_type != Periodic) { 
	int ibnd = static_cast<int>(bc_type); // either LO_DIRICHLET or LO_NEUMANN
	const Geometry& geom = bd.getGeom();

	// We are on the low side of the domain in coordinate direction n
	if (bx.smallEnd(n) == geom.Domain().smallEnd(n)) {
	  // Set the boundary conditions to live exactly on the faces of the domain
	  bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.0 );
	  
	  // Set the Dirichlet/Neumann boundary values 
	  bd.setValue(Orientation(n, Orientation::low) ,i, bc_value);
	  
	  // Define the type of boundary conditions 
	  bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,ibnd);
	}
	
	// We are on the high side of the domain in coordinate direction n
	if (bx.bigEnd(n) == geom.Domain().bigEnd(n)) {
	  // Set the boundary conditions to live exactly on the faces of the domain
	  bd.setBoundLoc(Orientation(n, Orientation::high) ,i,0.0 );
	  
	  // Set the Dirichlet/Neumann boundary values
	  bd.setValue(Orientation(n, Orientation::high) ,i, bc_value);
	  
	  // Define the type of boundary conditions 
	  bd.setBoundCond(Orientation(n, Orientation::high) ,i,comp,ibnd);
	}
      } 
    }
  }
}

void solve(MultiFab& soln, const MultiFab& anaSoln, 
	   Real a, Real b, MultiFab& alpha, MultiFab beta[], MultiFab& beta_cc,
	   MultiFab& rhs, const BoxArray& bs, const Geometry& geom,
	   solver_t solver)
{
  BL_PROFILE("solve()");
  std::string ss;

  soln.setVal(0.0);

  const Real run_strt = ParallelDescriptor::second();

  if (solver == BoxLib_C) {
    ss = "CPP";
    solve_with_Cpp(soln, a, b, alpha, beta, rhs, bs, geom);
  }
#ifdef USE_F90_SOLVERS
  else if (solver == BoxLib_F) {
    ss = "F90";
    solve_with_F90(soln, a, b, alpha, beta, rhs, bs, geom);
  }
#endif
#ifdef USEHYPRE
  else if (solver == Hypre) {
    ss = "Hyp";
    solve_with_hypre(soln, a, b, alpha, beta, rhs, bs, geom);
  }
#endif
#ifdef USEHPGMG
  else if (solver == HPGMG) {
    ss = "HPGMG";
    solve_with_HPGMG(soln, a, b, alpha, beta_cc, rhs, bs, geom, n_cell);
  }
#endif
  else {
    BoxLib::Error("Invalid solver");
  }

  Real run_time = ParallelDescriptor::second() - run_strt;

  ParallelDescriptor::ReduceRealMax(run_time, ParallelDescriptor::IOProcessorNumber());
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "   Run time      : " << run_time << std::endl;
  }

  if (plot_soln) {
    writePlotFile("SOLN-"+ss, soln, geom);
  }

  if (plot_err || comp_norm) {
    soln.minus(anaSoln, 0, Ncomp, 0); // soln contains errors now
    MultiFab& err = soln;

    if (plot_err) {
      writePlotFile("ERR-"+ss, soln, geom);
    }
    
    if (comp_norm) {
      Real twoNorm = err.norm2();
      Real maxNorm = err.norm0();
      
      err.setVal(1.0);
      Real vol = err.norm2();
      twoNorm /= vol;
      
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "   2 norm error  : " << twoNorm << std::endl;
	std::cout << "   max norm error: " << maxNorm << std::endl;
      }
    }
  }
}

void solve_with_Cpp(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab beta[], 
		    MultiFab& rhs, const BoxArray& bs, const Geometry& geom)
{
  BL_PROFILE("solve_with_Cpp()");
  BndryData bd(bs, 1, geom);
  set_boundary(bd, rhs);

  ABecLaplacian abec_operator(bd, dx);
  abec_operator.setScalars(a, b);
  abec_operator.setCoefficients(alpha, beta);

  MultiGrid mg(abec_operator);
  mg.setVerbose(verbose);
  mg.solve(soln, rhs, tolerance_rel, tolerance_abs);
} 

#ifdef USE_F90_SOLVERS
void solve_with_F90(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab beta[], 
		    MultiFab& rhs, const BoxArray& bs, const Geometry& geom)
{
  BL_PROFILE("solve_with_F90()");
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

  // Set the boundary conditions to live exactly on the faces of the domain --
  //   this is only relevant for Dirichlet and Neumann bc's
  for ( int n = 0; n < BL_SPACEDIM; ++n ) {
    xa[0][n] = 0.;
    xb[0][n] = 0.;
  }

  if (bc_type == Periodic) {
    // Define the type of boundary conditions to be periodic
    for ( int n = 0; n < BL_SPACEDIM; ++n ) {
      mg_bc[2*n + 0] = MGT_BC_PER;
      mg_bc[2*n + 1] = MGT_BC_PER;
    }
  }
  else if (bc_type == Neumann) {
    // Define the type of boundary conditions to be Neumann
    for ( int n = 0; n < BL_SPACEDIM; ++n ) {
      mg_bc[2*n + 0] = MGT_BC_NEU;
      mg_bc[2*n + 1] = MGT_BC_NEU;
    }
  }
  else if (bc_type == Dirichlet) {
    // Define the type of boundary conditions to be Dirichlet
    for ( int n = 0; n < BL_SPACEDIM; ++n ) {
      mg_bc[2*n + 0] = MGT_BC_DIR;
      mg_bc[2*n + 1] = MGT_BC_DIR;
    }
  }

  MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false, CC_CROSS_STENCIL);

  MultiFab* soln_p[1]; soln_p[0] = &soln;
  MultiFab*  rhs_p[1];  rhs_p[0] = &rhs;

  PArray<MultiFab> acoeffs(1, PArrayManage);
  acoeffs.set(0, new MultiFab(bs,Ncomp,0,Fab_allocate));
  acoeffs[0].copy(alpha);
  acoeffs[0].mult(a); 

  Array< PArray<MultiFab> > bcoeffs(BL_SPACEDIM);
  for (int n = 0; n < BL_SPACEDIM ; n++) {
    bcoeffs[n].resize(1,PArrayNoManage);
    bcoeffs[n].set(0, &beta[n]);
  }

  // The coefficients are set such that we will solve
  //  (a alpha - b del dot beta grad) soln = rhs
  //  written in the form 
  //  (acoeffs - b del dot bcoeffs grad) soln = rhs
  mgt_solver.set_visc_coefficients(acoeffs,bcoeffs,b,xa,xb);

  BCRec phys_bc;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    phys_bc.setLo(n,0);
    phys_bc.setHi(n,0);
  }

  MacBndry bndry(bs,1,geom);
  bndry.setBndryValues(soln,0,0,1,phys_bc);
  
  int always_use_bnorm = 0;
  Real final_resnorm;
  mgt_solver.solve(soln_p, rhs_p, bndry, tolerance_rel, tolerance_abs, always_use_bnorm, final_resnorm);
}
#endif

#ifdef USEHYPRE
void solve_with_hypre(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab beta[], 
		      MultiFab& rhs, const BoxArray& bs, const Geometry& geom)
{
  BL_PROFILE("solve_with_hypre()");
  BndryData bd(bs, 1, geom);
  set_boundary(bd, rhs);

  HypreABecLap hypreSolver(bs, geom);
  hypreSolver.setScalars(a, b);
  hypreSolver.setACoeffs(alpha);
  hypreSolver.setBCoeffs(beta);
  hypreSolver.setVerbose(verbose);
  hypreSolver.solve(soln, rhs, tolerance_rel, tolerance_abs, maxiter, bd);
}
#endif

#ifdef USEHPGMG
void solve_with_HPGMG(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab& beta_cc,
		      MultiFab& rhs, const BoxArray& bs, const Geometry& geom, int n_cell)
{
  int minCoarseDim;
  if (domain_boundary_condition == BC_PERIODIC)
  {
    minCoarseDim = 2; // avoid problems with black box calculation of D^{-1} for poisson with periodic BC's on a 1^3 grid
  }
  else
  {
    minCoarseDim = 1; // assumes you can drop order on the boundaries
  }

  level_type level_h;
  mg_type MG_h;
  int numVectors = 12;

  int my_rank = 0, num_ranks = 1;

#ifdef BL_USE_MPI
  MPI_Comm_size (MPI_COMM_WORLD, &num_ranks);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
#endif /* BL_USE_MPI */

  const double h0 = dx[0];
  // Create the geometric structure of the HPGMG grid using the RHS MultiFab as
  // a template. This doesn't copy any actual data.
  MultiFab::CreateHPGMGLevel(&level_h, rhs, n_cell, max_grid_size, my_rank, num_ranks, domain_boundary_condition, numVectors, h0);

  // Set up the coefficients for the linear operator L.
  MultiFab::SetupHPGMGCoefficients(a, b, alpha, beta_cc, &level_h);

  // Now that the HPGMG grid is built, populate it with RHS data.
  MultiFab::ConvertToHPGMGLevel(rhs, n_cell, max_grid_size, &level_h, VECTOR_F);

#ifdef USE_HELMHOLTZ
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Creating Helmholtz (a=" << a << ", b=" << b << ") test problem" << std::endl;;
  }
#else
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Creating Poisson (a=" << a << ", b=" << b << ") test problem" << std::endl;;
  }
#endif /* USE_HELMHOLTZ */

  if (level_h.boundary_condition.type == BC_PERIODIC)
  {
    double average_value_of_f = mean (&level_h, VECTOR_F);
    if (average_value_of_f != 0.0)
    {
      if (ParallelDescriptor::IOProcessor())
      {
        std::cerr << "WARNING: Periodic boundary conditions, but f does not sum to zero... mean(f)=" << average_value_of_f << std::endl;
      }
      //shift_vector(&level_h,VECTOR_F,VECTOR_F,-average_value_of_f);
    }
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rebuild_operator(&level_h,NULL,a,b);    // i.e. calculate Dinv and lambda_max
  MGBuild(&MG_h,&level_h,a,b,minCoarseDim); // build the Multigrid Hierarchy
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (ParallelDescriptor::IOProcessor())
      std::cout << std::endl << std::endl << "===== STARTING SOLVE =====" << std::endl << std::flush;

  MGResetTimers (&MG_h);
  zero_vector (MG_h.levels[0], VECTOR_U);
#ifdef USE_FCYCLES
  FMGSolve (&MG_h, 0, VECTOR_U, VECTOR_F, a, b, tolerance_abs, tolerance_rel);
#else
  MGSolve (&MG_h, 0, VECTOR_U, VECTOR_F, a, b, tolerance_abs, tolerance_rel);
#endif /* USE_FCYCLES */

  MGPrintTiming (&MG_h, 0);   // don't include the error check in the timing results
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (ParallelDescriptor::IOProcessor())
    std::cout << std::endl << std::endl << "===== Performing Richardson error analysis ==========================" << std::endl;
  // solve A^h u^h = f^h
  // solve A^2h u^2h = f^2h
  // solve A^4h u^4h = f^4h
  // error analysis...
  MGResetTimers(&MG_h);
  const double dtol = tolerance_abs;
  const double rtol = tolerance_rel;
  int l;for(l=0;l<3;l++){
    if(l>0)restriction(MG_h.levels[l],VECTOR_F,MG_h.levels[l-1],VECTOR_F,RESTRICT_CELL);
           zero_vector(MG_h.levels[l],VECTOR_U);
    #ifdef USE_FCYCLES
    FMGSolve(&MG_h,l,VECTOR_U,VECTOR_F,a,b,dtol,rtol);
    #else
     MGSolve(&MG_h,l,VECTOR_U,VECTOR_F,a,b,dtol,rtol);
    #endif
  }
  richardson_error(&MG_h,0,VECTOR_U);

  // Now convert solution from HPGMG back to rhs MultiFab.
  MultiFab::ConvertFromHPGMGLevel(soln, &level_h, VECTOR_U);

  const double norm_from_HPGMG = norm(&level_h, VECTOR_U);
  const double mean_from_HPGMG = mean(&level_h, VECTOR_U);
  const Real norm0 = soln.norm0();
  const Real norm2 = soln.norm2();
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "mean from HPGMG: " << mean_from_HPGMG << std::endl;
    std::cout << "norm from HPGMG: " << norm_from_HPGMG << std::endl;
    std::cout << "norm0 of RHS copied to MF: " << norm0 << std::endl;
    std::cout << "norm2 of RHS copied to MF: " << norm2 << std::endl;
  }

  // Write the MF to disk for comparison with the in-house solver
  if (plot_soln)
  {
    writePlotFile("SOLN-HPGMG", soln, geom);
  }

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  MGDestroy(&MG_h);
  destroy_level(&level_h);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}
#endif
