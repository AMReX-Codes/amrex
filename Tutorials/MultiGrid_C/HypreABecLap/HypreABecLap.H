#ifndef _HypreABecLap_H_
#define _HypreABecLap_H_

#include <iomanip>
#include <iostream>

#include <Tuple.H>
#include <Geometry.H>
#include <MultiFab.H>
#include <MacBndry.H>
#include <ParmParse.H>

#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"

class HypreABecLap {

public:

  HypreABecLap(const BoxArray& grids, const Geometry& geom,
	       int solver_flag = 0);
  ~HypreABecLap();

  void setScalars(Real sa, Real sb);
  void setACoeffs(const MultiFab& alpha);
  void setBCoeffs(const MultiFab beta[]);
  void setVerbose(int _verbose);
  void solve(MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
	     int max_iter, const BndryData& bndry);

protected:

  int solver_flag, verbose;

  int is_periodic[BL_SPACEDIM];
  Real dx[BL_SPACEDIM];

  MultiFab* acoefs;
  Tuple<MultiFab*, BL_SPACEDIM> bcoefs;

  Real scalar_a, scalar_b;

  int maxiter;
  Real reltol, abstol;

  HYPRE_StructGrid grid;
  HYPRE_StructMatrix A;
  HYPRE_StructVector b, x;
  HYPRE_StructSolver solver, precond;

  int pfmg_relax_type;
};
#endif

