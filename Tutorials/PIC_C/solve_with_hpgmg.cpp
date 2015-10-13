#ifdef USEHPGMG
#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#include <stencil_types.H>

#include "Particles.H"

void
Gravity::solve_with_HPGMG(MultiFab& rhs,
                          MultiFab& grad_phi,
                          const Geometry& geom,
                          Real tol,
                          Real abs_tol)
{
  const Geometry& geom = parent->Geom(level);
  const Real* dx = parent->Geom(level).CellSize();
  const Box& domain = geom.Domain();

  MultiFab soln(rhs.boxArray(),1,1);

  BndryData bd(grids[level], 1, geom);
  set_boundary(bd, rhs, dx);
  const int n_cell = domain.length(0);

  // Because HPGMG solves the MINUS Laplacian. Isn't this fun?
  rhs.mult(-1.0, 0, 1, 0);

  // We'll get the max grid (box) size from the soln MultiFab already provided
  // to us. Just pick the first valid Box we can find and measure its size.
  // HPGMG requires the Boxes to be cubes, so if they're not then a sanity
  // check in MultiFab::CreateHPGMGLevel() will catch it and quit.
  MFIter mfi(soln);
  while (!mfi.isValid()) ++mfi;
  const Box& bx = mfi.validbox();
  const int max_grid_size = bx.length(0);

  const int my_rank = ParallelDescriptor::MyProc();
  const int num_ranks = ParallelDescriptor::NProcs();

  Laplacian lap_operator(bd, dx[0]);
  const Real a = 0.0;
  const Real b = 1.0;

  const BoxArray& ba = rhs.boxArray();
  MultiFab alpha(ba, 1, 0, Fab_allocate);
  MultiFab beta_cc(ba, 1, 1, Fab_allocate);
  alpha.setVal(0.0);
  beta_cc.setVal(1.0);

  const Real tolerance_abs = 0.0;
  const Real tolerance_rel = 1.0e-10;

  const int domain_boundary_condition = BC_PERIODIC;

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

  MultiFab::CreateHPGMGLevel(&level_h, rhs, n_cell, max_grid_size, my_rank, num_ranks, domain_boundary_condition, numVectors, *dx);
  MultiFab::SetupHPGMGCoefficients(a, b, alpha, beta_cc, &level_h);
  MultiFab::ConvertToHPGMGLevel(rhs, n_cell, max_grid_size, &level_h, VECTOR_F);

  if (level_h.boundary_condition.type == BC_PERIODIC)
  {
    Real average_value_of_f = mean (&level_h, VECTOR_F);
    if (average_value_of_f != 0.0)
    {
      if (ParallelDescriptor::IOProcessor())
      {
        std::cerr << "WARNING: Periodic boundary conditions, but f does not sum to zero... mean(f)=" << average_value_of_f << std::endl;
      }
      shift_vector(&level_h,VECTOR_F,VECTOR_F,-average_value_of_f);
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

  if (show_timings)
    MGPrintTiming (&MG_h, 0);

  // Now convert solution from HPGMG back to rhs MultiFab.
  MultiFab::ConvertFromHPGMGLevel(soln, &level_h, VECTOR_U);

  MGDestroy(&MG_h); // destroys all but the finest grid
  destroy_level(&level_h); // destroys the finest grid

  soln.FillBoundary();
  geom.FillPeriodicBoundary(soln);

  lap_operator.compFlux(grad_phi[0],grad_phi[1],grad_phi[2],soln);
}
#endif
