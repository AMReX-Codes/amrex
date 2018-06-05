#ifndef AMREX_HypreABecLap_H_
#define AMREX_HypreABecLap_H_

#include <iomanip>
#include <iostream>

#include <AMReX_Tuple.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MacBndry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_LO_BCTYPES.H>

#include <HYPRE_struct_ls.h>

namespace amrex {

class HypreABecLap {

public:

  HypreABecLap(const BoxArray& grids,
               const DistributionMapping& dmap,
               const Geometry& geom_,
               MPI_Comm comm_);
  ~HypreABecLap();

  void setScalars(Real sa, Real sb);
  void setACoeffs(const MultiFab& alpha);
  void setBCoeffs(const MultiFab beta[]);
  void setBCoeffs(const std::array<const MultiFab*,BL_SPACEDIM>& beta);
  void setVerbose(int _verbose);
  void solve(MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
	     int max_iter, const BndryData& bndry);
  void solve(MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
	     int max_iter, LinOpBCType bc_type, Real bc_value);

protected:

  MPI_Comm comm;

  int solver_flag, verbose;

  int is_periodic[BL_SPACEDIM];
  Real dx[BL_SPACEDIM];

  amrex::Geometry geom;

  MultiFab* acoefs;
  Tuple<MultiFab*, BL_SPACEDIM> bcoefs;

  Real scalar_a, scalar_b;

  int maxiter;
  Real reltol, abstol;

  HYPRE_StructGrid grid;
  HYPRE_StructMatrix A;
  HYPRE_StructVector b, x;
  HYPRE_StructSolver solver, precond;

  int pfmg_rap_type, pfmg_relax_type, skip_relax;
  int num_pre_relax, num_post_relax;
  int print_level;

  // static utility functions follow:

  static int ispow2(int i)
  {
    return (i == 1) ? 1 : (((i <= 0) || (i & 1)) ? 0 : ispow2(i / 2));
  }
  
  static int* loV(const Box& b) {
    return (int*) b.loVect();
  }
  
  static int* hiV(const Box& b) {
    return (int*) b.hiVect();
  }

};

}

#endif

