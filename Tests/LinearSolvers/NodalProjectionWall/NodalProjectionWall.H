#ifndef NODAL_PROJECTION_WALL_H_
#define NODAL_PROJECTION_WALL_H_

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

// MLNodeLaplacian::setWallDivuCorrection for every compRHS and projection of the test (input
// wall_corr)
extern bool wall_corr;

// Burggraf flow: u = 8 f(x) g'(y), v = -8 f'(x) g(y), f = x^4 - 2x^3 + x^2, g = y^4 - y^2
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real bf  (amrex::Real x) { return x*x*x*x - amrex::Real(2)*x*x*x + x*x; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real bf1 (amrex::Real x) { return amrex::Real(4)*x*x*x - amrex::Real(6)*x*x + amrex::Real(2)*x; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real bf2 (amrex::Real x) { return amrex::Real(12)*x*x - amrex::Real(12)*x + amrex::Real(2); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real bf3 (amrex::Real x) { return amrex::Real(24)*x - amrex::Real(12); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real bg  (amrex::Real y) { return y*y*y*y - y*y; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real bg1 (amrex::Real y) { return amrex::Real(4)*y*y*y - amrex::Real(2)*y; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real bg2 (amrex::Real y) { return amrex::Real(12)*y*y - amrex::Real(2); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real bg3 (amrex::Real y) { return amrex::Real(24)*y; }

// The checks of compRHS without a solve (checks.cpp).  Prints them and returns the number of
// failed checks.
int run_checks (amrex::Vector<int> const& n_cell, int max_grid_size);

// The convergence study with projections and viscous operators (study.cpp)
void run_study (amrex::Vector<int> const& n_cell, int max_grid_size, amrex::Real rtol, int verbose);

#endif
