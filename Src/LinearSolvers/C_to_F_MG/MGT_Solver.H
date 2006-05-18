#ifndef _MGT_SOLVER_

#define _MGT_SOLVER_

#include <MultiFab.H>
#include <Geometry.H>
#include <BCRec.H>

#include <mg_cpp_f.h>

class MGT_Solver
{
public:
  MGT_Solver(int nlevel, const BoxArray& ba, const Box& domain, const BCRec& phys_bc, const double* dx,
	     const DistributionMapping& dmap, const Geometry& geom);
  void solve(MultiFab& uu, MultiFab& rh);
  ~MGT_Solver();
private:
  int m_mgt;
  const int m_dm;
  const int m_nlevel;
  const BoxArray m_ba;
  const Box m_domain;
  const DistributionMapping& m_dmap;
  const Geometry& m_geom;
  const BCRec& m_phys_rec;
  double m_dx[BL_SPACEDIM];
};
#endif
