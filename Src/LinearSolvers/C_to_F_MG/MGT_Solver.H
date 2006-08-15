#ifndef _MGT_SOLVER_
#define _MGT_SOLVER_

#include <vector>

#include <MultiFab.H>
#include <Geometry.H>
#include <BCRec.H>
#include <BndryData.H>

#include <mg_cpp_f.h>

class MGT_Solver
{
public:
      
  MGT_Solver(const std::vector<Geometry>& geom,
             const BCRec& phys_bc, 
	     const std::vector<BoxArray>& grids, 
	     const std::vector<DistributionMapping>& dmap,
	     bool nodal);
  static void initialize(bool nodal);

  void set_coefficients(const MultiFab* aa[], const MultiFab* bb[][BL_SPACEDIM],
                        const BndryData& bd);
  void solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
             const BndryData& bd);

  void set_nodal_coefficients(const MultiFab* aa[]);
  void nodal_project(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol);

  ~MGT_Solver();

  static int def_maxiter, def_maxiter_b, def_bottom_solver;
  static int def_nu_1, def_nu_2, def_nu_b, def_nu_f, def_gamma;
  static Real def_omega, def_bottom_solver_eps;
  static int def_verbose, def_cg_verbose, def_min_width, def_max_nlevel;
  static int def_cycle, def_smoother;
  static int def_usecg, def_cg_solver;
  static int stencil_type;
  
private:

  MultiFab* m_acoefs;
  MultiFab* m_bcoefs[BL_SPACEDIM];
  Real m_alpha;
  Real m_beta;
  int m_nlevel;

  const std::vector<DistributionMapping> m_dmap;
  const std::vector<BoxArray> m_grids;
  const bool m_nodal;

  static bool initialized;

};
#endif
