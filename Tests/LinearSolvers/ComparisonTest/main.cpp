// We solve (a alpha - b del dot beta grad) soln = rhs
// where a and b are scalars, alpha and beta are arrays

#include <Utility.H>
#include <ParallelDescriptor.H>
#include <ParmParse.H>
#include <PArray.H>

#include <LO_BCTYPES.H>

#ifdef USEHYPRE
#include <HypreABecLap.H>
#endif

#include <COEF_F.H>
#include <writePlotFile.H>

enum solver_t {BoxLib_C, BoxLib_F, Hypre, All};
enum bc_t {Periodic = 0,
	   Dirichlet = LO_DIRICHLET, 
	   Neumann = LO_NEUMANN};
solver_t solver_type = All;
bc_t     bc_type = Periodic;

void build_grids(std::vector<Geometry>& geom, 
		 std::vector<BoxArray>& grids);
void setup_coef(PArray<MultiFab> &exac, PArray<MultiFab> &alph, 
		PArray<MultiFab> &beta, PArray<MultiFab> &rhs, 
		const std::vector<Geometry>& geom, 
		const std::vector<BoxArray>& grids,
		Real a, Real b, Real sigma, Real w);
void solve_with_F90(PArray<MultiFab>& soln, Real a, Real b, 
		    const PArray<MultiFab>& alph, 
		    const PArray<MultiFab>& beta, 
		    PArray<MultiFab>& rhs, 
		    const std::vector<Geometry>& geom, 
		    const std::vector<BoxArray>& grids,
		    int ibnd);
#ifdef USEHYPRE
void solve_with_hypre(PArray<MultiFab>& soln, Real a, Real b, 
		      const PArray<MultiFab>& alph, 
		      const PArray<MultiFab>& beta, 
		      PArray<MultiFab>& rhs, 
		      const std::vector<Geometry>& geom, 
		      const std::vector<BoxArray>& grids,
		      int ibnd);
#endif
void compute_norm(const PArray<MultiFab>& soln, const PArray<MultiFab>& exac, 
		  const std::vector<Geometry>& geom, const std::vector<BoxArray>& grids,
		  int nsoln, int iCpp, int iF90, int iHyp);

int main(int argc, char* argv[])
{
  BoxLib::Initialize(argc,argv);

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
      BoxLib::Error("");
    }
  }

  int max_level = 0;
  pp.query("max_level", max_level);
  int nlevel = max_level+1;

  std::vector<Geometry> geom(nlevel);
  std::vector<BoxArray> grids(nlevel);

  build_grids(geom, grids);

  PArray<MultiFab> soln(nlevel, PArrayManage);
  PArray<MultiFab> soln1(nlevel, PArrayNoManage);
  PArray<MultiFab> exac(nlevel, PArrayManage);
  PArray<MultiFab> alph(nlevel, PArrayManage);
  PArray<MultiFab> beta(nlevel, PArrayManage);
  PArray<MultiFab> rhs(nlevel, PArrayManage);

  int nsoln=-1, iF90=-1, iCpp=-1, iHyp=-1;
  switch (solver_type) 
    {
    case BoxLib_C:
      nsoln = 1;
      iCpp = 0;
      break;
    case BoxLib_F:
      nsoln = 1;
      iF90 = 0;
      break;
    case Hypre:
      nsoln = 1;
      iHyp = 0;
      break;
    case All:
#ifdef USEHYPRE
      nsoln = 3;
      iCpp = 0;
      iF90 = 1;
      iHyp = 2;
#else
      nsoln = 2;
      iCpp = 0;
      iF90 = 1;
#endif
      break;
    }

  for (int ilev=0; ilev < nlevel; ilev++) {
    soln.set(ilev, new MultiFab(grids[ilev], nsoln, 1));
    if (nsoln == 1) {
      soln1.set(ilev, &(soln.get(ilev)));
    }
    else {
      soln1.set(ilev, new MultiFab(grids[ilev], 1, 1));
    }
    exac.set(ilev, new MultiFab(grids[ilev], 1, 0));
    alph.set(ilev, new MultiFab(grids[ilev], 1, 0));
    beta.set(ilev, new MultiFab(grids[ilev], 1, 1)); // one ghost cell
    rhs.set (ilev, new MultiFab(grids[ilev], 1, 0));
  }

  Real a, b, sigma, w;
  pp.get("a",  a);
  pp.get("b",  b);
  pp.get("sigma", sigma);
  pp.get("w", w);  

  setup_coef(exac, alph, beta, rhs, geom, grids, a, b, sigma, w);

  int ibnd = static_cast<int>(bc_type);

  if (solver_type == BoxLib_C || solver_type == All) {
    for (int ilev=0; ilev < nlevel; ilev++) {
      soln1[ilev].setVal(0.0);
    }    

    //    solve_with_Cpp(soln1, a, b, alph, beta, rhs, geom, grids, nlevel, ibnd);

    if (nsoln > 1) { // soln1 doesn't point to the same multifabs as soln
      for (int ilev=0; ilev < nlevel; ilev++) {
	MultiFab::Copy(soln[ilev], soln1[ilev], 0, iCpp, 1, 1);
      }        
    }
  }

  if (solver_type == BoxLib_F || solver_type == All) {
    for (int ilev=0; ilev < nlevel; ilev++) {
      soln1[ilev].setVal(0.0);
    }    

    solve_with_F90(soln1, a, b, alph, beta, rhs, geom, grids, ibnd);

    if (nsoln > 1) { // soln1 doesn't point to the same multifabs as soln
      for (int ilev=0; ilev < nlevel; ilev++) {
	MultiFab::Copy(soln[ilev], soln1[ilev], 0, iF90, 1, 1);
      }        
    }
  }

#ifdef USEHYPRE
  if (solver_type == Hypre || solver_type == All) {
    for (int ilev=0; ilev < nlevel; ilev++) {
      soln1[ilev].setVal(0.0);
    }    

    solve_with_hypre(soln1, a, b, alph, beta, rhs, geom, grids, ibnd);

    if (nsoln > 1) { // soln1 doesn't point to the same multifabs as soln
      for (int ilev=0; ilev < nlevel; ilev++) {
	MultiFab::Copy(soln[ilev], soln1[ilev], 0, iHyp, 1, 1);
      }        
    }
  }
#endif

  if (nsoln > 1) { // soln1 doesn't point to the same multifabs as soln,
                   // and soln1 is defined with PArrayNoManage
    for (int ilev=0; ilev < nlevel; ilev++) {
      delete &(soln1[ilev]);
    }
  }

  int write_plot = 0;
  pp.query("write_plot", write_plot);
  if (write_plot) {
    writePlotFile("plot", soln, exac, alph, beta, rhs, geom, grids, nsoln, iCpp, iF90, iHyp);
  }

  int comp_norm = 1;
  pp.query("comp_norm", comp_norm);
  if (comp_norm) {
    compute_norm(soln, exac, geom, grids, nsoln, iCpp, iF90, iHyp);
  }

  BoxLib::Finalize();
}

void build_grids(std::vector<Geometry>& geom, std::vector<BoxArray>& grids)
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

void setup_coef(PArray<MultiFab> &exac, PArray<MultiFab> &alph, 
		PArray<MultiFab> &beta, PArray<MultiFab> &rhs, 
		const std::vector<Geometry>& geom, 
		const std::vector<BoxArray>& grids,
		Real a, Real b, Real sigma, Real w)
{
  int ibnd = static_cast<int>(bc_type); 
  int nlevel = geom.size();

  for (int ilev=0; ilev < nlevel; ilev++) {
    const Geometry& geo = geom[ilev];
    const Real* dx = geo.CellSize();  

    for (MFIter mfi(alph[ilev]); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      const Box& bx = grids[ilev][i];
    
      FORT_SET_COEF(exac[ilev][mfi].dataPtr(), alph[ilev][mfi].dataPtr(),
		    beta[ilev][mfi].dataPtr(), rhs[ilev][mfi].dataPtr(),
		    bx.loVect(), bx.hiVect(), geo.ProbLo(), geo.ProbHi(),
		    dx, a, b, sigma, w, ibnd);
    }
  }
}
