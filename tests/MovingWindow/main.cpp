
#include <BoxLib.H>
#include <ParmParse.H>
#include <Array.H>
#include <MultiFab.H>
#include <PlotFileUtil.H>

#include <WarpXConst.H>

#include "movingWindow_F.H"

void testMovingWindow() {
  
  int n_steps, n_cell, max_grid_size, is_periodic[BL_SPACEDIM];
  Array<Real> moving_window_x(BL_SPACEDIM, 0.0);
  Array<Real> moving_window_v(BL_SPACEDIM, 0.0);
  Array<Real> prob_lo(BL_SPACEDIM);
  Array<Real> prob_hi(BL_SPACEDIM);

  Real dt = 1.2e-15;

  // parse inputs
  {
    ParmParse pp;
    pp.get("n_steps", n_steps);
    pp.get("n_cell", n_cell);
    pp.get("max_grid_size", max_grid_size);
    pp.getarr("prob_lo", prob_lo, 0, BL_SPACEDIM);
    pp.getarr("prob_hi", prob_hi, 0, BL_SPACEDIM);
    pp.getarr("prob_lo", moving_window_x, 0, SPACEDIM);
    pp.getarr("moving_window_v", moving_window_v, 0, BL_SPACEDIM);
  }

  // multiply the moving window velocity by the speed of light
  for (int dir = 0; dir < BL_SPACEDIM; dir++) {
    moving_window_v[dir] *= PhysConst::c;
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
  int Nghost = 1;

  MultiFab E(ba, 1, Nghost);
  E.setVal(0.0);

  for ( MFIter mfi(E); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.validbox();

    init_E(E[mfi].dataPtr(),
	   bx.loVect(), bx.hiVect(), &Nghost,
	   geom.CellSize(), geom.ProbLo(), geom.ProbHi());
  }

  // take some "time steps"
  for (int step = 0; step < n_steps; step++) {

    MultiFab outputE(ba, 1, 0);
    MultiFab::Copy(outputE, E, 0, 0, 1, 0);
    const std::string& plotname = BoxLib::Concatenate("plt", step, 5);
    Array<std::string> varnames{"E"};
    BoxLib::WriteSingleLevelPlotfile(plotname, outputE, varnames, geom, 0.0, 0);

    for (int dir=0; dir<BL_SPACEDIM; dir++) {
      moving_window_x[dir] += moving_window_v[dir] * dt;
    }

    // shift E one cell here.
    for ( MFIter mfi(E); mfi.isValid(); ++mfi )
      {
	const Box& bx = mfi.validbox();
	
	shift_E(E[mfi].dataPtr(),
		bx.loVect(), bx.hiVect(), &Nghost,
		geom.CellSize(), geom.ProbLo(), geom.ProbHi());
      }

    E.FillBoundary(geom.periodicity());

  }   
}

int main(int argc, char* argv[])
{
  BoxLib::Initialize(argc,argv);

  testMovingWindow();

  BoxLib::Finalize();
}
