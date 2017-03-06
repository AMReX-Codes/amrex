// We solve (a alpha - b del dot beta grad) soln = rhs
// where a and b are scalars, alpha and beta are arrays

#include <fstream>
#include <iomanip>

#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_BndryData.H>
#include <AMReX_MultiGrid.H>
#include <AMReX_CGSolver.H>
#include <AMReX_Laplacian.H>
#include <AMReX_ABecLaplacian.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MacBndry.H>

#include <COEF_F.H>
#include <RHS_F.H>

using namespace amrex;

Real tolerance_rel = 1.e-8;
Real tolerance_abs = 0.0;
int  maxiter       = 1; 
int  comp_norm     = 1;
int  verbose       = 0;

Real dx[BL_SPACEDIM];

enum solver_t {BoxLib_C, BoxLib_F, All};
enum bc_t {Periodic = 0,
	   Dirichlet = LO_DIRICHLET, 
	   Neumann = LO_NEUMANN};
solver_t solver_type;
bc_t     bc_type;

int Ncomp = 1;

const int maxGrid(64);

void compute_analyticSolution(MultiFab& anaSoln);
void setup_coeffs(BoxArray& bs, MultiFab& alpha, MultiFab beta[], const Geometry& geom);
void setup_rhs(MultiFab& rhs, const Geometry& geom, Real a, Real b);
void set_boundary(BndryData& bd, const MultiFab& rhs, int comp);
void solve(MultiFab& soln, const MultiFab& anaSoln, 
	   Real a, Real b, MultiFab& alpha, MultiFab beta[], 
	   MultiFab& rhs, const BoxArray& bs, const Geometry& geom,
	   solver_t solver);


int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);

  BL_PROFILE_VAR("main()", pmain);

  std::cout << std::setprecision(15);

  solver_type = BoxLib_C;
  bc_type = Periodic;

  Real     a = 0.0;
  Real     b = 1.0;

  // ---- First use the number of processors to decide how many grids you have.
  // ---- We arbitrarily decide to have one grid per MPI process in a uniform
  // ---- cubic domain, so we require that the number of processors be N^3.
  // ---- This requirement is somewhat arbitrary, but convenient for now.

  int nprocs = ParallelDescriptor::NProcs();

  // N is the cube root of the number of processors
  int N(0);
  for(int i(1); i*i*i <= nprocs; ++i) {
    if(i*i*i == nprocs) {
      N = i;
    }
  }

  if(N == 0) {  // not a cube
    if(ParallelDescriptor::IOProcessor()) {
      std::cerr << "**** Error:  nprocs = " << nprocs << " is not currently supported." << std::endl;
    }
    amrex::Error("We require that the number of processors be a perfect cube");
  }
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "N = " << N << std::endl;
  }


  // ---- make a box, then a boxarray with maxSize
  int domain_hi = (N*maxGrid) - 1;
  Box domain(IntVect(0,0,0), IntVect(domain_hi,domain_hi,domain_hi));
  BoxArray bs(domain);
  bs.maxSize(maxGrid);

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
    for (int n = 0; n < BL_SPACEDIM; n++) is_per[n] = 0;
  } 
  else {
    for (int n = 0; n < BL_SPACEDIM; n++) is_per[n] = 1;
  }
 
  // This defines a Geometry object which is useful for writing the plotfiles
  Geometry geom(domain,&real_box,coord,is_per);

  for ( int n=0; n<BL_SPACEDIM; n++ ) {
    dx[n] = ( geom.ProbHi(n) - geom.ProbLo(n) )/domain.length(n);
  }

  if (ParallelDescriptor::IOProcessor()) {
     std::cout << "Domain size     : " << N << std::endl;
     std::cout << "Max_grid_size   : " << maxGrid << std::endl;
     std::cout << "Number of grids : " << bs.size() << std::endl;
  }

  DistributionMapping dm{bs};

  // Allocate and define the right hand side.
  MultiFab rhs(bs, dm, Ncomp, 0); 
  setup_rhs(rhs, geom, a, b);

  MultiFab alpha(bs, dm, Ncomp, 0);
  MultiFab beta[BL_SPACEDIM];
  for ( int n=0; n<BL_SPACEDIM; ++n ) {
    BoxArray bx(bs);
    beta[n].define(bx.surroundingNodes(n), dm, Ncomp, 1);
  }

  setup_coeffs(bs, alpha, beta, geom);

  MultiFab anaSoln;
  if (comp_norm) {
      anaSoln.define(bs, dm, Ncomp, 0);
      compute_analyticSolution(anaSoln);
  }

  // Allocate the solution array 
  // Set the number of ghost cells in the solution array.
  MultiFab soln(bs, dm, Ncomp, 1);

  solve(soln, anaSoln, a, b, alpha, beta, rhs, bs, geom, BoxLib_C);

  BL_PROFILE_VAR_STOP(pmain);

  amrex::Finalize();
}

void compute_analyticSolution(MultiFab& anaSoln)
{
  BL_PROFILE("compute_analyticSolution");

  int ibnd = static_cast<int>(bc_type); 

  for (MFIter mfi(anaSoln); mfi.isValid(); ++mfi) {
    const int* alo = anaSoln[mfi].loVect();
    const int* ahi = anaSoln[mfi].hiVect();
    const Box& bx = mfi.validbox();

    FORT_COMP_ASOL(anaSoln[mfi].dataPtr(), ARLIM(alo), ARLIM(ahi),
		   bx.loVect(),bx.hiVect(),dx, ibnd);
  }
}

void setup_coeffs(BoxArray& bs, MultiFab& alpha, MultiFab beta[], const Geometry& geom)
{
  BL_PROFILE("setup_coeffs");

  Real sigma = 1.0;
  Real     w = 1.0;

  MultiFab cc_coef(bs,alpha.DistributionMap(), Ncomp,1); // cell-centered beta

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
}

void setup_rhs(MultiFab& rhs, const Geometry& geom, Real a, Real b)
{
  BL_PROFILE("setup_rhs");
  Real sigma = 1.0;
  Real     w = 1.0;

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
}

void set_boundary(BndryData& bd, const MultiFab& rhs, int comp=0)
{
  BL_PROFILE("set_boundary");
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
	   Real a, Real b, MultiFab& alpha, MultiFab beta[], 
	   MultiFab& rhs, const BoxArray& bs, const Geometry& geom,
	   solver_t solver)
{
  BL_PROFILE("solve");
  soln.setVal(0.0);

  const Real run_strt = ParallelDescriptor::second();

  BndryData bd(bs, soln.DistributionMap(), 1, geom);
  set_boundary(bd, rhs);

  ABecLaplacian abec_operator(bd, dx);
  abec_operator.setScalars(a, b);
  abec_operator.setCoefficients(alpha, beta);

  MultiGrid mg(abec_operator);
  mg.setMaxIter(maxiter);
  mg.setVerbose(verbose);
  mg.setFixedIter(1);
  mg.solve(soln, rhs, tolerance_rel, tolerance_abs);

  Real run_time = ParallelDescriptor::second() - run_strt;

  ParallelDescriptor::ReduceRealMax(run_time, ParallelDescriptor::IOProcessorNumber());

  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Run time        : " << run_time << std::endl;
  }
}
