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
  MGT_Solver(const BndryData& bd, 
	     const double* dx, 
	     const std::vector<BoxArray>& grids, 
	     const std::vector<DistributionMapping>& dmap,
	     bool nodal);
  void solve(MultiFab* uu[], MultiFab* rh[]);
  ~MGT_Solver();
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
};
#endif
