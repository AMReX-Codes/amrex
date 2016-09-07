#ifndef _MGT_SOLVER_
#define _MGT_SOLVER_

#include <vector>
#include <MultiFab.H>
#include <Geometry.H>
#include <BCRec.H>
#include <BndryData.H>
#include <PArray.H>

#include <mg_cpp_f.h>

class MGT_Solver
{
public:

    MGT_Solver(const std::vector<Geometry>& geom, 
               int* bc, 
               const std::vector<BoxArray>& grids,
               const std::vector<DistributionMapping>& dmap,
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
  static void set_cfa2 (const MultiFab& aa, int lev);
  static void set_cfb (const MultiFab& bb, Real beta, int lev, int dir);
  static void set_cfbn (const MultiFab& bb, Real beta, int lev, int dir);

  static void set_rh (const MultiFab& mf, int lev);
  static void set_uu (const MultiFab& mf, int lev);

  static void get_uu (MultiFab& mf, int lev, int ng);
  static void get_res (MultiFab& mf, int lev);
  static void get_gp (MultiFab& mf, int lev, int dir, Real dx);

  //
  // (alpha * aa - beta * (del dot bb grad)) phi = RHS
  //
  void set_abeclap_coeffs (Real alpha, // aa is constant one
			   Real beta,
			   const Array<PArray<MultiFab> >& bb,
			   const Array< Array<Real> >& xa,
			   const Array< Array<Real> >& xb);
  void set_abeclap_coeffs (const PArray<MultiFab>& aa, // alpha is one
			   Real beta,
			   const Array<PArray<MultiFab> >& bb,
			   const Array< Array<Real> >& xa,
			   const Array< Array<Real> >& xb);
  void set_abeclap_coeffs (Real alpha,
			   const PArray<MultiFab>& aa,
			   Real beta,
			   const Array<PArray<MultiFab> >& bb,
			   const Array< Array<Real> >& xa,
			   const Array< Array<Real> >& xb);

  void set_mac_coefficients(const Array<PArray<MultiFab> >& bb,
                            const Array< Array<Real> >& xa,
                            const Array< Array<Real> >& xb);

  void set_gravity_coefficients(const Array< PArray<MultiFab> >& bb,
                                const Array< Array<Real> >& xa,
                                const Array< Array<Real> >& xb);

  void set_const_gravity_coeffs(const Array< Array<Real> >& xa,
                                const Array< Array<Real> >& xb);

  void set_porous_coefficients(const PArray<MultiFab>& a1, 
			       const PArray<MultiFab>& a2, 
			       const Array<PArray<MultiFab> >& bb, 
			       const Real& beta, 
			       const Array< Array<Real> >& xa, 
			       const Array< Array<Real> >& xb,
			       int nc_opt = 0);

  void set_maxorder(const int max_order);

  void solve(MultiFab* uu[], MultiFab* rh[], const BndryData& bd,
	     Real tol, Real abs_tol, int always_use_bnorm, 
	     Real& final_resnorm, int need_grad_phi=0);

  void applyop(MultiFab* uu[], MultiFab* res[], const BndryData& bd);

  void compute_residual(MultiFab* uu[], MultiFab* rh[], MultiFab* res[], const BndryData& bd);

  void get_fluxes(int lev, PArray<MultiFab>& flux, const Real* dx);

  void set_nodal_coefficients(const MultiFab* aa[]);
  void set_nodal_const_coefficients(Real val);
  void nodal_project(MultiFab* uu[], MultiFab* vel[], MultiFab* rhcc[], const PArray<MultiFab>& rhnd, 
		     const Real& tol, const Real& abs_tol, 
		     int* lo_inflow=0, int* hi_inflow=0);

  void fill_sync_resid(MultiFab* sync_resid_crse, const MultiFab& msk, const MultiFab& vold,
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

   void Build(const std::vector<Geometry>& geom, 
	      int* bc, int  stencil_type,
	      const std::vector<DistributionMapping>& dmap,
	      int nc,
	      int ncomp);

  int verbose;
  int m_nlevel;
  std::vector<BoxArray> m_grids;
  bool m_nodal;
  bool have_rhcc;

  static bool initialized;

};
#endif
