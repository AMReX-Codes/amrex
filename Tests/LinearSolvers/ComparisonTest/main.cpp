// We solve (a alpha - b del dot beta grad) soln = rhs
// where a and b are scalars, alpha and beta are arrays

#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <AMReX_LO_BCTYPES.H>

#ifdef USEHYPRE
#include <HypreABecLap.H>
#endif

#include <COEF_F.H>
#include <writePlotFile.H>

using namespace amrex;

enum solver_t {BoxLib_C, BoxLib_F, Hypre, All};
enum bc_t {Periodic = 0,
	   Dirichlet = LO_DIRICHLET, 
	   Neumann = LO_NEUMANN};
solver_t solver_type = All;
bc_t     bc_type = Periodic;

void build_grids(Vector<Geometry>& geom, 
		 Vector<BoxArray>& grids);
void setup_coef(const Vector<MultiFab*> &exac,
		const Vector<MultiFab*> &alph, 
		const Vector<MultiFab*> &beta,
		const Vector<MultiFab*> &rhs, 
		const Vector<Geometry>& geom, 
		const Vector<BoxArray>& grids,
		Real a, Real b, Real sigma, Real w);
#ifdef USEHYPRE
void solve_with_hypre(const Vector<MultiFab*>& soln, Real a, Real b, 
		      const Vector<MultiFab*>& alph, 
		      const Vector<MultiFab*>& beta, 
		      const Vector<MultiFab*>& rhs, 
		      const Vector<Geometry>& geom, 
		      const Vector<BoxArray>& grids,
		      int ibnd);
#endif
void compute_norm(const Vector<MultiFab*>& soln,
		  const Vector<MultiFab*>& exac, 
		  const Vector<Geometry>& geom, const Vector<BoxArray>& grids,
		  int nsoln, int iCpp, int iHyp);

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);

  {

  ParmParse pp;

  {
    std::string solver_type_s;
    pp.get("solver_type",solver_type_s);
    if (solver_type_s == "BoxLib_C") {
      solver_type = BoxLib_C;
    }
    else if (solver_type_s == "BoxLib_F") {
      solver_type = BoxLib_F;      
    }
    else if (solver_type_s == "Hypre") {
#ifdef USEHYPRE
      solver_type = Hypre;
#else
      amrex::Error("Set USE_HYPRE=TRUE in GNUmakefile");
#endif
    }
    else if (solver_type_s == "All") {
      solver_type = All;
    }  
    else {
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Don't know this solver type: " << solver_type << std::endl;
      }
      amrex::Error("");
    }
  }

  {
    std::string bc_type_s;
    pp.get("bc_type",bc_type_s);
    if (bc_type_s == "Dirichlet") {
      bc_type = Dirichlet;
    }
    else if (bc_type_s == "Neumann") {
      bc_type = Neumann;
    }
    else if (bc_type_s == "Periodic") {
      bc_type = Periodic;
    }
    else {
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Don't know this boundary type: " << bc_type << std::endl;
      }
      amrex::Error("");
    }
  }

  int max_level = 0;
  pp.query("max_level", max_level);
  int nlevel = max_level+1;

  Vector<Geometry> geom(nlevel);
  Vector<BoxArray> grids(nlevel);

  build_grids(geom, grids);

  Vector<std::unique_ptr<MultiFab> >  soln(nlevel);
  Vector<std::unique_ptr<MultiFab> > soln1(nlevel);
  Vector<std::unique_ptr<MultiFab> >  exac(nlevel);
  Vector<std::unique_ptr<MultiFab> >  alph(nlevel);
  Vector<std::unique_ptr<MultiFab> >  beta(nlevel);
  Vector<std::unique_ptr<MultiFab> >   rhs(nlevel);

  int nsoln=-1, iCpp=-1, iHyp=-1;
  switch (solver_type) 
    {
    case BoxLib_C:
      nsoln = 1;
      iCpp = 0;
      break;
    case Hypre:
      nsoln = 1;
      iHyp = 0;
      break;
    case All:
#ifdef USEHYPRE
      nsoln = 3;
      iCpp = 0;
      iHyp = 2;
#else
      nsoln = 2;
      iCpp = 0;
#endif
      break;
    }

  for (int ilev=0; ilev < nlevel; ilev++) {
    DistributionMapping dm {grids[ilev]};
    soln [ilev].reset(new MultiFab(grids[ilev], dm, nsoln, 1));
    soln1[ilev].reset(new MultiFab(grids[ilev], dm, 1, 1));
    exac [ilev].reset(new MultiFab(grids[ilev], dm, 1, 0));
    alph [ilev].reset(new MultiFab(grids[ilev], dm, 1, 0));
    beta [ilev].reset(new MultiFab(grids[ilev], dm, 1, 1)); // one ghost cell
    rhs  [ilev].reset(new MultiFab(grids[ilev], dm, 1, 0));
  }

  auto psoln  = amrex::GetVecOfPtrs(soln);
  auto psoln1 = amrex::GetVecOfPtrs(soln1);
  auto pexac  = amrex::GetVecOfPtrs(exac);
  auto palph  = amrex::GetVecOfPtrs(alph);
  auto pbeta  = amrex::GetVecOfPtrs(beta);
  auto prhs   = amrex::GetVecOfPtrs(rhs);

  Real a, b, sigma, w;
  pp.get("a",  a);
  pp.get("b",  b);
  pp.get("sigma", sigma);
  pp.get("w", w);  

  setup_coef(pexac, palph, pbeta, prhs, geom, grids, a, b, sigma, w);

  int ibnd = static_cast<int>(bc_type);

  if (solver_type == BoxLib_C || solver_type == All) {
    for (int ilev=0; ilev < nlevel; ilev++) {
	soln1[ilev]->setVal(0.0);
    }    

    //    solve_with_Cpp(soln1, a, b, alph, beta, rhs, geom, grids, nlevel, ibnd);

    for (int ilev=0; ilev < nlevel; ilev++) {
	MultiFab::Copy(*soln[ilev], *soln1[ilev], 0, iCpp, 1, 1);
    }
  }

#ifdef USEHYPRE
  if (solver_type == Hypre || solver_type == All) {
    for (int ilev=0; ilev < nlevel; ilev++) {
	soln1[ilev]->setVal(0.0);
    }    

    solve_with_hypre(psoln1, a, b, palph, pbeta, prhs, geom, grids, ibnd);

    for (int ilev=0; ilev < nlevel; ilev++) {
	MultiFab::Copy(*soln[ilev], *soln1[ilev], 0, iHyp, 1, 1);
    }
  }
#endif

  int write_plot = 0;
  pp.query("write_plot", write_plot);
  if (write_plot) {
      writePlotFile("plot", psoln, pexac, palph, pbeta, prhs, geom, grids, nsoln, iCpp, iHyp);
  }

  int comp_norm = 1;
  pp.query("comp_norm", comp_norm);
  if (comp_norm) {
      compute_norm(psoln, pexac, geom, grids, nsoln, iCpp, iHyp);
  }

  }

  amrex::Finalize();
}

void build_grids(Vector<Geometry>& geom, Vector<BoxArray>& grids)
{
  ParmParse pp;
  
  // build the gridse level
  int n_cell;
  int max_grid_size;

  pp.get("n_cell",n_cell);
  pp.get("max_grid_size",max_grid_size);

  // Define a single box covering the domain
#if (BL_SPACEDIM == 1)
  IntVect dom0_lo(0);
  IntVect dom0_hi(n_cell-1);
#endif
#if (BL_SPACEDIM == 2)
  IntVect dom0_lo(0,0);
  IntVect dom0_hi(n_cell-1,n_cell-1);
#endif
#if (BL_SPACEDIM == 3)
  IntVect dom0_lo(0,0,0);
  IntVect dom0_hi(n_cell-1,n_cell-1,n_cell-1);
#endif

  Box dom0(dom0_lo,dom0_hi);

  BoxArray ba0(dom0);

  grids[0] = ba0; 
  grids[0].maxSize(max_grid_size);

  int nlevel=grids.size();
  for (int ilev=1; ilev<nlevel; ilev++) {
    ba0.grow(-n_cell/4);
    ba0.refine(2);

    grids[ilev] = ba0;
    grids[ilev].maxSize(max_grid_size);    
  }

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
 
  geom[0].define(dom0, &real_box, coord, is_per);
  for (int ilev=1; ilev<nlevel; ilev++) {
    dom0.refine(2);
    geom[ilev].define(dom0, &real_box, coord, is_per);
  }
}

void setup_coef(const Vector<MultiFab*> &exac,
		const Vector<MultiFab*> &alph, 
		const Vector<MultiFab*> &beta,
		const Vector<MultiFab*> &rhs, 
		const Vector<Geometry>& geom, 
		const Vector<BoxArray>& grids,
		Real a, Real b, Real sigma, Real w)
{
  int ibnd = static_cast<int>(bc_type); 
  int nlevel = geom.size();

  for (int ilev=0; ilev < nlevel; ilev++) {
    const Geometry& geo = geom[ilev];
    const Real* dx = geo.CellSize();  

    for (MFIter mfi(*alph[ilev]); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      const Box& bx = grids[ilev][i];
    
      FORT_SET_COEF((*exac[ilev])[mfi].dataPtr(), (*alph[ilev])[mfi].dataPtr(),
		    (*beta[ilev])[mfi].dataPtr(), (*rhs[ilev])[mfi].dataPtr(),
		    bx.loVect(), bx.hiVect(), geo.ProbLo(), geo.ProbHi(),
		    dx, a, b, sigma, w, ibnd);
    }
  }
}
