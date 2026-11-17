#ifndef NODAL_PROJECTION_WALL_H_
#define NODAL_PROJECTION_WALL_H_

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

// MLNodeLaplacian::setBoundaryDivuCorrection for every compRHS of the test (input wall_corr)
extern bool wall_corr;

// The checks of compRHS without a solve (checks.cpp).  Prints them and returns the number of
// failed checks.
int run_checks (amrex::Vector<int> const& n_cell, int max_grid_size);

#endif
