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
  enum 
    {
      // IPAR address: match to mg_cpp.f90
      MGT_NU1 = 0,
      MGT_NU2 = 1,
      MGT_GAMMA = 2,
      MGT_MAX_ITER = 3,
      MGT_MAX_NLEVEL = 4,
      MGT_MIN_WIDTH  = 5,
      MGT_VERBOSE = 6,
      MGT_BOTTOM_SOLVER = 7,
      MGT_BOTTOM_MAX_ITER = 8,
      MGT_SMOOTHER = 9,
      MGT_CYCLE = 10,
      // RPAR ADDRESS: match to mg_cpp.f90
      MGT_EPS  = 0,
      MGT_OMEGA = 1,
      MGT_BOTTOM_SOLVER_EPS = 2
    };
      
  MGT_Solver(const BndryData& bd, 
	     const BCRec& phys_bc, 
	     const double* dx, 
	     const std::vector<BoxArray>& grids, 
	     const std::vector<DistributionMapping>& dmap,
	     const std::vector<int>& ipar, const std::vector<double>& rpar,
	     bool nodal);
  void set_coefficients(const MultiFab* aa[], const MultiFab* bb[][BL_SPACEDIM]);
  void solve(MultiFab* uu[], MultiFab* rh[]);
  ~MGT_Solver();

  static std::vector<int> ipar_defaults();
  static std::vector<double> rpar_defaults();
  
private:
  int m_mgt;

  MultiFab* m_acoefs;
  MultiFab* m_bcoefs[BL_SPACEDIM];
  Real m_alpha;
  Real m_beta;
  double m_dx[BL_SPACEDIM];
  int m_nlevel;

  const BndryData& m_bd;
  const std::vector<DistributionMapping> m_dmap;
  const std::vector<BoxArray> m_grids;
  const bool m_nodal;
  std::vector<int> m_ipar;
  std::vector<double> m_rpar;
};
#endif
