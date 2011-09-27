// We solve (a alpha - b del dot beta grad) soln = rhs
// where a and b are scalars, alpha and beta are arrays

#include <Utility.H>
#include <ParallelDescriptor.H>
#include <ParmParse.H>
#include <PArray.H>

#include <LO_BCTYPES.H>

#include "COEF_F.H"
#include "writePlotFile.H"

enum solver_t {BoxLib_C, BoxLib_F, All};
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
void solve_with_F90(PArray<MultiFab>& soln, int iF90, Real a, Real b, 
		    const PArray<MultiFab>& alph, 
		    const PArray<MultiFab>& beta, 
		    PArray<MultiFab>& rhs, 
		    const std::vector<Geometry>& geom, 
		    const std::vector<BoxArray>& grids,
		    int ibnd);
void compute_norm(const PArray<MultiFab>& soln, const PArray<MultiFab>& exac, 
		  const std::vector<Geometry>& geom, const std::vector<BoxArray>& grids,
		  int nsoln, int iF90, int iCpp);

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
  PArray<MultiFab> exac(nlevel, PArrayManage);
  PArray<MultiFab> alph(nlevel, PArrayManage);
  PArray<MultiFab> beta(nlevel, PArrayManage);
  PArray<MultiFab> rhs(nlevel, PArrayManage);

  int nsoln=-1, iF90=-1, iCpp=-1;
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
    case All:
      nsoln = 2;
      iCpp = 0;
      iF90 = 1;
      break;
    }

  for (int ilev=0; ilev < nlevel; ilev++) {
    soln.set(ilev, new MultiFab(grids[ilev], nsoln, 1));
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

  for (int ilev=0; ilev < nlevel; ilev++) {
    soln[ilev].setVal(0.0);
  }    

  if (solver_type == BoxLib_C || solver_type == All) {
    //    solve_with_Cpp(soln, iCpp, a, b, alph, beta, rhs, geom, grids, nlevel, ibnd);
  }
  else if (solver_type == BoxLib_F || solver_type == All) {
    int ibnd = static_cast<int>(bc_type);
    solve_with_F90(soln, iF90, a, b, alph, beta, rhs, geom, grids, ibnd);
  }

  int write_plot = 0;
  pp.query("write_plot", write_plot);
  if (write_plot) {
    writePlotFile("plot", soln, exac, alph, beta, rhs, geom, grids, nsoln, iF90, iCpp);
  }

  int comp_norm = 1;
  pp.query("comp_norm", comp_norm);
  if (comp_norm) {
    compute_norm(soln, exac, geom, grids, nsoln, iF90, iCpp);
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
  IntVect dom0_lo(0,0,0);
  IntVect dom0_hi(n_cell-1,n_cell-1,n_cell-1);
  Box dom0(dom0_lo,dom0_hi);

  BoxArray ba0(dom0);

  grids[0].define(ba0); 
  grids[0].maxSize(max_grid_size);

  int nlevel=grids.size();
  for (int ilev=1; ilev<nlevel; ilev++) {
    ba0.grow(-n_cell/4);
    ba0.refine(2);

    grids[ilev].define(ba0);
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
  int* is_per = new int[BL_SPACEDIM];
  
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
    
      FORT_SET_COEF(exac[ilev][i].dataPtr(), alph[ilev][i].dataPtr(),
		    beta[ilev][i].dataPtr(), rhs[ilev][i].dataPtr(),
		    bx.loVect(), bx.hiVect(), geo.ProbLo(), geo.ProbHi(),
		    dx, a, b, sigma, w, ibnd);
    }
  }
}
