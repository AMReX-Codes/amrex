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
               bool verbose = false);

    void Build(const std::vector<Geometry>& geom, 
               int* bc, int  stencil_type,
               const std::vector<DistributionMapping>& dmap,
               int nc,
               int ncomp,
               bool verbose);

  static void initialize(bool nodal);
  static void Finalize();
  static void FlushFortranOutput();

  void set_mac_coefficients(const MultiFab* aa[], const MultiFab* bb[][BL_SPACEDIM],
                            Array< Array<Real> >& xa,
                            Array< Array<Real> >& xb);

  void set_gravity_coefficients(Array< PArray<MultiFab> >& coeffs,
                                Array< Array<Real> >& xa,
                                Array< Array<Real> >& xb,
                                int is_constant);

  void set_const_gravity_coeffs(Array< Array<Real> >& xa,
                                Array< Array<Real> >& xb);

  // If index_order == 0, bb[dim][level]
  // else,                bb[level][dim] 
  // The second way is consistent with set_gravity_coefficients.
  void set_visc_coefficients(PArray<MultiFab>& aa, 
			     Array<PArray<MultiFab> >& bb, 
			     const Real& beta, 
			     Array< Array<Real> >& xa,
			     Array< Array<Real> >& xb, 
			     int index_order=0);

  void set_visc_coefficients(const MultiFab* aa[], 
                             const MultiFab* bb[][BL_SPACEDIM],
                             const Real& b, 
                             Array< Array<Real> >& xa,
                             Array< Array<Real> >& xb);

  void set_visc_coefficients(MultiFab* aa[], 
                             MultiFab* bb[][BL_SPACEDIM],
                             const Real& b, 
                             Array< Array<Real> >& xa,
                             Array< Array<Real> >& xb);

  void set_porous_coefficients(PArray<MultiFab>& a1, 
			       PArray<MultiFab>& a2, 
			       Array<PArray<MultiFab> >& bb, 
			       const Real& beta, 
			       Array< Array<Real> >& xa, 
			       Array< Array<Real> >& xb,
			       int nc_opt = 0);

  void set_porous_coefficients(const MultiFab* a1[], const MultiFab* a2[], 
                               const MultiFab* bb[][BL_SPACEDIM],
                               const Real& b, 
                               Array< Array<Real> >& xa,
                               Array< Array<Real> >& xb,
			       int nc_opt = 0);

  void set_porous_coefficients(MultiFab* a1[], const MultiFab* a2[], 
			       MultiFab* bb[][BL_SPACEDIM], 
			       const Real& beta, 
			       Array< Array<Real> >& xa,
			       Array< Array<Real> >& xb,
			       int nc_opt = 0);

  void set_maxorder(const int max_order);

  // *******************************************************************************
  // These two take a status flag and pass the status flag to mgt_solve
  // *******************************************************************************

  void solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
             const BndryData& bd, int need_grad_phi, Real& final_resnorm,int& status);

  void solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
	     int need_grad_phi, Real& final_resnorm,int& status);

  // *******************************************************************************
  // These four do *not* take a status flag but do pass an internal status flag to mgt_solve
  // *******************************************************************************

  void solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
             const BndryData& bd, Real& final_resnorm);

  void solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
             const BndryData& bd, int need_grad_phi, Real& final_resnorm);

  void solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
	     const BndryData bd[], int need_grad_phi, Real& final_resnorm);

  void solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
	     int need_grad_phi, Real& final_resnorm);

  // *******************************************************************************
  // End of solve options
  // *******************************************************************************

  void applyop(MultiFab* uu[], MultiFab* res[], const BndryData& bd);

  void applyop(MultiFab* uu[], MultiFab* res[], const BndryData bd[]);

  void applybc(MultiFab* uu[], const BndryData& bd);

  void compute_residual(MultiFab* uu[], MultiFab* rh[], MultiFab* res[], const BndryData& bd);

  void get_fluxes(int lev, PArray<MultiFab>& flux, const Real* dx);

  void get_fluxes(int lev, 
		  MultiFab* const* flux, 
		  const Real* dx);

  void set_nodal_coefficients(const MultiFab* aa[]);
  void nodal_project(MultiFab* uu[], MultiFab* vel[], MultiFab* rhcc[], MultiFab* rhcrse[], 
		     const Real& tol, const Real& abs_tol, 
		     int* lo_inflow=0, int* hi_inflow=0);

  void fill_sync_resid(MultiFab* sync_resid_crse, const MultiFab& msk, const MultiFab& vold,
		       int isCoarse);

  ~MGT_Solver();

  static int def_maxiter, def_maxiter_b, def_bottom_solver;
  static int def_nu_1, def_nu_2, def_nu_b, def_nu_f, def_gamma;
  static Real def_omega, def_bottom_solver_eps;
  static Real def_max_L0_growth;
  static int def_verbose, def_cg_verbose, def_min_width, def_max_nlevel;
  static int def_cycle, def_smoother;
  static int def_usecg, def_cg_solver;
  
private:

  MultiFab* m_acoefs;
  MultiFab* m_bcoefs[BL_SPACEDIM];
  Real m_alpha;
  Real m_beta;
  int m_nlevel;

  std::vector<BoxArray> m_grids;
  bool m_nodal;
  bool have_rhcc;

  static bool initialized;

};
#endif
