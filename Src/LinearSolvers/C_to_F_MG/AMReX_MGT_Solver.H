#ifndef AMREX_MGT_SOLVER_
#define AMREX_MGT_SOLVER_

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>
#include <AMReX_BndryData.H>

#include <mg_cpp_f.h>

namespace amrex {

class MGT_Solver
{
public:

    MGT_Solver(const Vector<Geometry>& geom, 
               int* bc, 
               const Vector<BoxArray>& grids,
               const Vector<DistributionMapping>& dmap,
               bool nodal,
               int stencil_type,
               bool _have_rhcc = false,
               int nc = 0,
               int ncomp = 1,
               int verbose = 0);
    
  static void initialize(bool nodal);
  static void Finalize();
  static void FlushFortranOutput();

  void set_cfa_const (Real alpha, int lev);
  static void set_cfa (const MultiFab& aa, int lev);
  static void set_cfaa (const MultiFab& aa, Real alpha, int lev);
  static void set_cfb (const MultiFab& bb, Real beta, int lev, int dir);

  static void set_rh (const MultiFab& mf, int lev);
  static void set_uu (const MultiFab& mf, int lev);

  static void get_uu (MultiFab& mf, int lev, int ng);
  static void get_res (MultiFab& mf, int lev);
  static void get_gp (MultiFab& mf, int lev, int dir, Real dx);

  static void set_rh_nodal (const MultiFab& mf, int lev);
  static void set_uu_nodal (const MultiFab& mf, int lev);
  static void get_uu_nodal (MultiFab& mf, int lev, int ng);

  //
  // (alpha * aa - beta * (del dot bb grad)) phi = RHS
  //
  void set_abeclap_coeffs (Real alpha, // aa is constant one
			   Real beta,
			   const Vector<Vector<MultiFab*> >& bb,
			   const Vector<Vector<Real> >& xa,
			   const Vector<Vector<Real> >& xb);
  void set_abeclap_coeffs (const Vector<MultiFab*>& aa, // alpha is one
			   Real beta,
			   const Vector<Vector<MultiFab*> >& bb,
			   const Vector<Vector<Real> >& xa,
			   const Vector<Vector<Real> >& xb);
  void set_abeclap_coeffs (Real alpha,
			   const Vector<MultiFab*>& aa,
			   Real beta,
			   const Vector<Vector<MultiFab*> >& bb,
			   const Vector<Vector<Real> >& xa,
			   const Vector<Vector<Real> >& xb);

  void set_mac_coefficients(const Vector<Vector<MultiFab*> >& bb,
                            const Vector<Vector<Real> >& xa,
                            const Vector<Vector<Real> >& xb);

  void set_gravity_coefficients(const Vector<Vector<MultiFab*> >& bb,
                                const Vector<Vector<Real> >& xa,
                                const Vector<Vector<Real> >& xb);

  void set_const_gravity_coeffs(const Vector< Vector<Real> >& xa,
                                const Vector< Vector<Real> >& xb);

  void set_maxorder(const int max_order);

  void solve(const Vector<MultiFab*>& uu, const Vector<MultiFab*>& rh, const BndryData& bd,
	     Real tol, Real abs_tol, int always_use_bnorm, 
	     Real& final_resnorm, int need_grad_phi=0);

  void solve_nodal(const Vector<MultiFab*>& uu, const Vector<MultiFab*>& rh,
                   Real tol, Real abs_tol);

  void applyop(const Vector<MultiFab*>& uu, const Vector<MultiFab*>& res, const BndryData& bd);

  void compute_residual(const Vector<MultiFab*>& uu, const Vector<MultiFab*>& rh,
			const Vector<MultiFab*>& res, const BndryData& bd);

  void get_fluxes(int lev, const Vector<MultiFab*>& flux, const Real* dx);

  void set_nodal_coefficients(const Vector<MultiFab*>& aa);
  void set_nodal_const_coefficients(Real val);
  void nodal_project(const Vector<MultiFab*>& uu, const Vector<MultiFab*>& vel, 
		     const Vector<MultiFab*>& rhcc, const Vector<MultiFab*>& rhnd, 
		     const Real& tol, const Real& abs_tol, 
		     int* lo_inflow=0, int* hi_inflow=0);

  void fill_sync_resid(MultiFab& sync_resid_crse, const MultiFab& msk, const MultiFab& vold,
		       int isCoarse);

  ~MGT_Solver();

  static int def_maxiter, def_maxiter_b, def_bottom_solver;
  static int def_nu_1, def_nu_2, def_nu_b, def_nu_f;
  static int def_verbose, def_cg_verbose;
  static int def_min_width, def_max_nlevel;
  static int def_cycle, def_smoother;
  static int def_usecg, def_cg_solver;
  static Real def_bottom_solver_eps, def_max_L0_growth;
  
private:

   void Build(const Vector<Geometry>& geom, 
	      int* bc, int  stencil_type,
	      const Vector<DistributionMapping>& dmap,
	      int nc,
	      int ncomp);

  int verbose;
  int m_nlevel;
  Vector<BoxArray> m_grids;
  Vector<DistributionMapping> m_dmap;
  bool m_nodal;
  bool have_rhcc;

  static bool initialized;

};

}

#endif
