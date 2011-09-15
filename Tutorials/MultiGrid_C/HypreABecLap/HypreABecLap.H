#ifndef _HypreABecLap_H_
#define _HypreABecLap_H_

#include <iomanip>
#include <iostream>

#include <Tuple.H>
#include <Geometry.H>
#include <MultiFab.H>
#include <MacBndry.H>
#include <ParmParse.H>

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

  int pfmg_relax_type, skip_relax;
  int num_pre_relax, num_post_relax;
  int print_level;

  // static utility functions follow:

  static int ispow2(int i)
  {
    return (i == 1) ? 1 : (((i <= 0) || (i & 1)) ? 0 : ispow2(i / 2));
  }
  
  static int* loV(const Box& b) {
#if (BL_SPACEDIM == 1)
    vl[0] = b.smallEnd(0);
    return vl;
#else
    return (int*) b.loVect();
#endif
  }
  
  static int* hiV(const Box& b) {
#if (BL_SPACEDIM == 1)
    vh[0] = b.bigEnd(0);
    return vh;
#else
    return (int*) b.hiVect();
#endif
  }

};
#endif

