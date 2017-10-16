
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <WarpXConst.H>
#include <WarpX.H>

#include "movingWindow_F.H"

using namespace amrex;

void testMovingWindow() {
  
  int n_steps, n_cell, max_grid_size, is_periodic[BL_SPACEDIM];
  int do_moving_window, moving_window_dir;
  Real moving_window_x = 0.0;
  Real moving_window_v = 0.0;
  Vector<Real> prob_lo(BL_SPACEDIM);
  Vector<Real> prob_hi(BL_SPACEDIM);

  Real dt = 1.2e-15;

  // parse inputs
  {
    ParmParse pp;
    pp.get("n_steps", n_steps);
    pp.get("n_cell", n_cell);
    pp.get("max_grid_size", max_grid_size);
    pp.get("do_moving_window", do_moving_window);
    std::string s;
    pp.get("moving_window_dir", s);
    if (s == "x" || s == "X") {
	moving_window_dir = 0;
    }
#if (BL_SPACEDIM == 3)
    else if (s == "y" || s == "Y") {
	moving_window_dir = 1;
    }
#endif
    else if (s == "z" || s == "Z") {
	moving_window_dir = BL_SPACEDIM-1;
    }
    else {
	const std::string msg = "Unknown moving_window_dir: "+s;
	amrex::Abort(msg.c_str()); 
    }
    pp.get("moving_window_v", moving_window_v);
    pp.getarr("prob_lo", prob_lo, 0, BL_SPACEDIM);
    pp.getarr("prob_hi", prob_hi, 0, BL_SPACEDIM);
  }

  moving_window_x = prob_lo[moving_window_dir];
  moving_window_v *= PhysConst::c;
  
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

    for (int i = 0; i < BL_SPACEDIM; i++) is_periodic[i] = 1;
    is_periodic[moving_window_dir] = 0;

    geom.define(domain, &real_box,coord, is_periodic);
  }

  const Real* dx = geom.CellSize();
  Real time = 0.0;
  int Nghost = 1;

  DistributionMapping dm {ba};
  MultiFab E(ba, dm, 1, Nghost);
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

    MultiFab outputE(ba, dm, 1, 0);
    MultiFab::Copy(outputE, E, 0, 0, 1, 0);
    const std::string& plotname = amrex::Concatenate("plt", step, 5);
    Vector<std::string> varnames{"E"};
    amrex::WriteSingleLevelPlotfile(plotname, outputE, varnames, geom, 0.0, 0);

    // update the window and figure out how much to shift
    int num_shift;
    int dir = moving_window_dir;

    Real new_lo[BL_SPACEDIM];
    Real new_hi[BL_SPACEDIM];    
    const Real* current_lo = geom.ProbLo();
    const Real* current_hi = geom.ProbHi();
    const Real* dx = geom.CellSize();

    moving_window_x += moving_window_v * dt;
    num_shift = (moving_window_x - current_lo[dir]) / dx[dir];

    for (int i=0; i<BL_SPACEDIM; i++) {
      new_lo[i] = current_lo[i];
      new_hi[i] = current_hi[i];
    }
    new_lo[dir] = current_lo[dir] + num_shift * dx[dir];
    new_hi[dir] = current_hi[dir] + num_shift * dx[dir];
    RealBox new_box(new_lo, new_hi);
    geom.ProbDomain(new_box);

    // shift the E Multifab
    WarpX::shiftMF(E, geom, num_shift, dir, IntVect::TheZeroVector());

  }
}

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);

  testMovingWindow();

  amrex::Finalize();
}
