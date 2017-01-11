
#include <random>

#include <BoxLib.H>
#include <ParmParse.H>
#include <Array.H>
#include <MultiFab.H>
#include <PlotFileUtil.H>

#include <WarpXConst.H>

void testMovingWindow() {
  
  int n_cell, max_grid_size, is_periodic[BL_SPACEDIM];
  Array<Real> moving_window_x(BL_SPACEDIM, 0.0);
  Array<Real> moving_window_v(BL_SPACEDIM, 0.0);
  Array<Real> prob_lo(BL_SPACEDIM);
  Array<Real> prob_hi(BL_SPACEDIM);

  // parse inputs
  {
    ParmParse pp;
    pp.get("n_cell", n_cell);
    pp.get("max_grid_size", max_grid_size);
    pp.getarr("prob_lo", prob_lo, 0, BL_SPACEDIM);
    pp.getarr("prob_hi", prob_hi, 0, BL_SPACEDIM);
    pp.getarr("moving_window_v", moving_window_v, 0, BL_SPACEDIM);
  }

  BoxArray ba;
  Geometry geom;
  {

    IntVect dom_lo(IntVect(D_DECL(0,0,0)));
    IntVect dom_hi(IntVect(D_DECL(n_cell-1, n_cell-1, n_cell-1)));
    Box domain(dom_lo, dom_hi);
    ba.define(domain);
    ba.maxSize(max_grid_size);

    RealBox real_box;
    real_box.setLo(prob_lo);
    real_box.setHi(prob_hi);

    int coord = 0;  // cartesian

    is_periodic[0] = 0;
    is_periodic[1] = 1;
    is_periodic[2] = 1;

    geom.define(domain, &real_box,coord, is_periodic);
  }

  const Real* dx = geom.CellSize();
  Real time = 0.0;
  
  MultiFab phi(ba, 1, 0);
  phi.setVal(0.0);

  std::string plotname{"plt00000"};
  Array<std::string> varnames{"phi"};
  BoxLib::WriteSingleLevelPlotfile(plotname, phi, varnames, geom, 0.0, 0);
}

int main(int argc, char* argv[])
{
  BoxLib::Initialize(argc,argv);

  testMovingWindow();

  BoxLib::Finalize();
}
