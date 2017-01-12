
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
    is_periodic[1] = 0;
    is_periodic[2] = 0;

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

    int num_shift[BL_SPACEDIM];
    Real new_lo[BL_SPACEDIM];
    Real new_hi[BL_SPACEDIM];
    
    const Real* current_lo = geom.ProbLo();
    const Real* current_hi = geom.ProbHi();
    const Real* dx = geom.CellSize();

    for (int dir=0; dir<BL_SPACEDIM; dir++) {
      moving_window_x[dir] += moving_window_v[dir] * dt;
      num_shift[dir] = (moving_window_x[dir] - current_lo[dir]) / dx[dir];
      new_lo[dir] = current_lo[dir] + num_shift[dir] * dx[dir];
      new_hi[dir] = current_hi[dir] + num_shift[dir] * dx[dir];
    }

    int max_shift = num_shift[0];
    for (int dir=1; dir<BL_SPACEDIM; dir++) {
      max_shift = std::max(max_shift, num_shift[dir]);
    }

    if (max_shift == 0) continue;
    if (max_shift <= Nghost) {

      if (num_shift[0] > 0) {
	for ( MFIter mfi(E); mfi.isValid(); ++mfi ) {
	  const Box& bx = mfi.validbox();
	  shift_x(E[mfi].dataPtr(),
		  bx.loVect(), bx.hiVect(), &Nghost, &num_shift[0]);
	}
      }

      if (num_shift[1] > 0) {
	for ( MFIter mfi(E); mfi.isValid(); ++mfi ) {
	  const Box& bx = mfi.validbox();
	  shift_y(E[mfi].dataPtr(),
		  bx.loVect(), bx.hiVect(), &Nghost, &num_shift[1]);
	}
      }

      if (num_shift[2] > 0) {
	for ( MFIter mfi(E); mfi.isValid(); ++mfi ) {
	  const Box& bx = mfi.validbox();
	  shift_z(E[mfi].dataPtr(),
		  bx.loVect(), bx.hiVect(), &Nghost, &num_shift[2]);
	}
      }
      E.FillBoundary(geom.periodicity());
    }

    if (max_shift > Nghost) {

      MultiFab tmpE(ba, 1, max_shift);
      MultiFab::Copy(tmpE, E, 0, 0, 1, 0);
      tmpE.FillBoundary(geom.periodicity());
      
      if (num_shift[0] > 0) {
	for ( MFIter mfi(tmpE); mfi.isValid(); ++mfi ) {
	  const Box& bx = mfi.validbox();
	  shift_x(tmpE[mfi].dataPtr(),
		  bx.loVect(), bx.hiVect(), &max_shift, &num_shift[0]);
	}
      }

      if (num_shift[1] > 0) {
	for ( MFIter mfi(tmpE); mfi.isValid(); ++mfi ) {
	  const Box& bx = mfi.validbox();
	  shift_y(tmpE[mfi].dataPtr(),
		  bx.loVect(), bx.hiVect(), &max_shift, &num_shift[1]);
	}
      }

      if (num_shift[2] > 0) {
	for ( MFIter mfi(tmpE); mfi.isValid(); ++mfi ) {
	  const Box& bx = mfi.validbox();
	  shift_z(tmpE[mfi].dataPtr(),
		  bx.loVect(), bx.hiVect(), &max_shift, &num_shift[2]);
	}
      }
      //      tmpE.FillBoundary(geom.periodicity());
      MultiFab::Copy(E, tmpE, 0, 0, 1, Nghost);
      E.FillBoundary(geom.periodicity());
    }

    RealBox new_box(new_lo, new_hi);
    geom.ProbDomain(new_box);

  }
}

int main(int argc, char* argv[])
{
  BoxLib::Initialize(argc,argv);

  testMovingWindow();

  BoxLib::Finalize();
}
