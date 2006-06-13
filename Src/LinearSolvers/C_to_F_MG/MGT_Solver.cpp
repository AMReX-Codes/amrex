#include <MGT_Solver.H>

typedef void (*mgt_get)(const int* mgt, const int* lev, const int* n, double* uu, 
			const int* plo, const int* phi, 
			const int* lo, const int* hi);
typedef void (*mgt_set)(const int* mgt, const int* lev, const int* n, const double* uu, 
			const int* plo, const int* phi, 
			const int* lo, const int* hi);
#if BL_SPACEDIM == 1
mgt_get mgt_get_uu = mgt_get_uu_1d;
mgt_set mgt_set_uu = mgt_set_uu_1d;
mgt_get mgt_get_rh = mgt_get_rh_1d;
mgt_set mgt_set_rh = mgt_set_rh_1d;
#elif BL_SPACEDIM == 2
mgt_get mgt_get_uu = mgt_get_uu_2d;
mgt_set mgt_set_uu = mgt_set_uu_2d;
mgt_get mgt_get_rh = mgt_get_rh_2d;
mgt_set mgt_set_rh = mgt_set_rh_2d;
#elif BL_SPACEDIM == 3
mgt_get mgt_get_uu = mgt_get_uu_3d;
mgt_set mgt_set_uu = mgt_set_uu_3d;
mgt_get mgt_get_rh = mgt_get_rh_3d;
mgt_set mgt_set_rh = mgt_set_rh_3d;
#endif
  
MGT_Solver::MGT_Solver(const BndryData& bd, const double* dx, 
		       const std::vector<BoxArray>& grids,
		       const std::vector<DistributionMapping>& dmap,
		       bool nodal)
  :
  m_bd(bd), m_dmap(dmap), m_grids(grids), m_nodal(nodal)
{
  BL_ASSERT( m_grids.size() == dmap.size() );
  m_nlevel = grids.size();
  int dm = BL_SPACEDIM;
  int i_nodal = (m_nodal)?1:0;
  mgt_alloc(&m_mgt, &dm, &m_nlevel, &i_nodal);
  for ( int i = 0; i < BL_SPACEDIM; ++i )
    {
      m_dx[i] = dx[i];
    }
  int nb = grids[0].size();
  std::vector<int> lo(nb*dm);
  std::vector<int> hi(nb*dm);
  for ( int i = 0; i < nb; ++i )
    {
      for ( int j = 0; j < dm; ++j )
	{
	  lo[i + j*nb] = m_grids[0][i].smallEnd(j);
	  hi[i + j*nb] = m_grids[0][i].bigEnd(j);
	}
    }
  int pm[dm];
  int bc[dm*2];
  const Geometry& geom = bd.getGeom();
  for ( int i = 0; i < dm; ++i ) 
    {
      pm[i] = geom.isPeriodic(i)? 1 : 0;
      if ( pm[i] )
	{
	  bc[i*2 + 0] = 0;
	  bc[i*2 + 1] = 0;
	}
      else
	{
	  // bc[i*2 + 0] = phys_bc.lo(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
	  // bc[i*2 + 1] = phys_bc.hi(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
	}
    }
  int lev = 0;
  const Array<int>& pmap = dmap[0].ProcessorMap();
  Box domain = bd.getDomain();
  mgt_set_level(&m_mgt, &lev, &nb, &dm, &lo[0], &hi[0], 
		domain.loVect(), domain.hiVect(), &bc[0], pm, &pmap[0]);
  mgt_finalize(&m_mgt);
}

#if 0
MGT_Solver::MGT_Solver(int nlevel, const BoxArray& ba, const Box& domain, const BCRec& phys_bc, const double* dx,
		       const DistributionMapping& dmap, const Geometry& geom)
  : m_dm(BL_SPACEDIM),
    m_nlevel(nlevel),
    m_ba(ba),
    m_domain(domain),
    m_dmap(dmap),
    m_geom(geom),
    m_phys_rec(phys_bc)
{
  mgt_alloc(&m_mgt, &m_dm, &m_nlevel);
  int nb = m_ba.size();
  std::vector<int> lo(nb*m_dm);
  std::vector<int> hi(nb*m_dm);
  for ( int i = 0; i < nb; ++i )
    {
      for ( int j = 0; j < m_dm; ++j )
	{
	  lo[i + j*nb] = m_ba[i].smallEnd(j);
	  hi[i + j*nb] = m_ba[i].bigEnd(j);
	}
    }
  int bc[m_dm*2];
  int pm[m_dm];
  for ( int i = 0; i < m_dm; ++i ) 
    {
      pm[i] = geom.isPeriodic(i)? 1 : 0;
      if ( pm[i] )
	{
	  bc[i*2 + 0] = 0;
	  bc[i*2 + 1] = 0;
	}
      else
	{
	  bc[i*2 + 0] = phys_bc.lo(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
	  bc[i*2 + 1] = phys_bc.hi(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
	}
    }
  int lev = 1;
  const Array<int>& pmap = dmap.ProcessorMap();
  mgt_set_level(&m_mgt, &lev, &nb, &m_dm, &lo[0], &hi[0], 
		m_domain.loVect(), m_domain.hiVect(), &bc[0], pm, &pmap[0]);
  const int verbose = 4;
  mgt_set_verbose(&m_mgt, &verbose);
  for ( int i = 0; i < m_dm; ++i )
    {
      m_dx[i] = dx[i];
    }
}
#endif

void 
MGT_Solver::solve(MultiFab* uu[], MultiFab* rh[])
{
  for ( int i = 0; i < m_nlevel; ++i )
    {
      int lev = i + 1;
      for (MFIter umfi(*(uu[i])); umfi.isValid(); ++umfi)
	{
	  const FArrayBox& rhs = (*(rh[i]))[umfi];
	  const FArrayBox& sol = (*(uu[i]))[umfi];
	  const Real* rd = rhs.dataPtr();
	  const Real* sd = sol.dataPtr();
	  int n = umfi.index() + 1; // Fortran Index Correction
	  const int* lo = umfi.validbox().loVect();
	  const int* hi = umfi.validbox().hiVect();
	  const int* plo = rhs.box().loVect();
	  const int* phi = rhs.box().hiVect();
	  mgt_set_rh(&m_mgt, &lev, &n, rd, plo, phi, lo, hi);
	  mgt_set_uu(&m_mgt, &lev, &n, sd, plo, phi, lo, hi);
	}
    }
  mgt_finalize(&m_mgt);
  double xa[BL_SPACEDIM], xb[BL_SPACEDIM];
  double pxa[BL_SPACEDIM], pxb[BL_SPACEDIM];
  for ( int i = 0; i < BL_SPACEDIM; ++i ) 
    {
      xa[i] = xb[i] = m_dx[i]/2;
      pxa[i] = pxb[i] = 0;
    }
  mgt_finalize_stencil(&m_mgt, xa, xb, pxa, pxb);
  mgt_solve_cc(&m_mgt);
  for ( int i = 0; i < m_nlevel; ++i )
    {
      int lev = i + 1;
      for (MFIter umfi(*(uu[i])); umfi.isValid(); ++umfi)
	{
	  FArrayBox& sol = (*(uu[i]))[umfi];
	  Real* sd = sol.dataPtr();
	  int n = umfi.index() + 1; // Fortran Index Correction
	  const int* lo = umfi.validbox().loVect();
	  const int* hi = umfi.validbox().hiVect();
	  const int* plo = sol.box().loVect();
	  const int* phi = sol.box().hiVect();
	  mgt_get_uu(&m_mgt, &lev, &n, sd, plo, phi, lo, hi);
	}
    }
}

MGT_Solver::~MGT_Solver()
{
  mgt_dealloc(&m_mgt);
}
