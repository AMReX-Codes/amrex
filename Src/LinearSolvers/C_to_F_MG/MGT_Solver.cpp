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
		       const std::vector<int>& ipar, const std::vector<double>& rpar,
		       bool nodal)
  :
  m_bd(bd), m_dmap(dmap), m_grids(grids), m_nodal(nodal), m_ipar(ipar), m_rpar(rpar)
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
	  bc[i*2 + 0] = MGT_BC_DIR;	// FIXME: Hardware DIRICHET
	  bc[i*2 + 1] = MGT_BC_DIR;
	  // bc[i*2 + 0] = phys_bc.lo(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
	  // bc[i*2 + 1] = phys_bc.hi(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
	}
    }
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      const Array<int>& pmap = dmap[lev].ProcessorMap();
      Box domain = bd.getDomain();
      mgt_set_level(&m_mgt, &lev, &nb, &dm, &lo[0], &hi[0], 
		    domain.loVect(), domain.hiVect(), &bc[0], pm, &pmap[0]);
    }
  int n_ipar = m_ipar.size();
  int n_rpar = m_rpar.size();
  double* rpar_p = n_rpar?&m_rpar[0]:0;
  int* ipar_p = n_ipar?&m_ipar[0]:0;
  mgt_finalize(&m_mgt, &n_ipar, ipar_p, &n_rpar, rpar_p);
}

void
MGT_Solver::set_coefficients(MultiFab* aa[], MultiFab* bb[][BL_SPACEDIM])
{
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      mgt_init_coeffs_lev(&m_mgt, &lev);
      double xa[BL_SPACEDIM], xb[BL_SPACEDIM];
      double pxa[BL_SPACEDIM], pxb[BL_SPACEDIM];
      for ( int i = 0; i < BL_SPACEDIM; ++i ) 
	{
	  xa[i] = xb[i] = m_dx[i]/2;
	  pxa[i] = pxb[i] = 0;
	}
      for (MFIter amfi(*(aa[lev])); amfi.isValid(); ++amfi)
	{
	  FArrayBox* a = &((*(aa[lev]))[amfi]);
	  FArrayBox* b[BL_SPACEDIM];
	  for ( int i = 0; i < BL_SPACEDIM; ++i )
	    {
	      b[i] = &((*(bb[lev][i]))[amfi]);
	    }
	}
      mgt_finalize_stencil_lev(&m_mgt, &lev, xa, xb, pxa, pxb);
    }
  mgt_finalize_stencil(&m_mgt);
}

void 
MGT_Solver::solve(MultiFab* uu[], MultiFab* rh[])
{
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      for (MFIter umfi(*(uu[lev])); umfi.isValid(); ++umfi)
	{
	  const FArrayBox& rhs = (*(rh[lev]))[umfi];
	  const FArrayBox& sol = (*(uu[lev]))[umfi];
	  const Real* rd = rhs.dataPtr();
	  const Real* sd = sol.dataPtr();
	  int n = umfi.index();
	  const int* lo = umfi.validbox().loVect();
	  const int* hi = umfi.validbox().hiVect();
	  const int* plo = rhs.box().loVect();
	  const int* phi = rhs.box().hiVect();
	  mgt_set_rh(&m_mgt, &lev, &n, rd, plo, phi, lo, hi);
	  mgt_set_uu(&m_mgt, &lev, &n, sd, plo, phi, lo, hi);
	}
    }
  mgt_solve(&m_mgt);
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      for (MFIter umfi(*(uu[lev])); umfi.isValid(); ++umfi)
	{
	  FArrayBox& sol = (*(uu[lev]))[umfi];
	  Real* sd = sol.dataPtr();
	  int n = umfi.index();
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

std::vector<int>
MGT_Solver::ipar_defaults()
{
  std::vector<int> ip(11);
  int n = ip.size();
  mgt_get_ipar_defaults(&ip[0], &n);
  return ip;
}

std::vector<double>
MGT_Solver::rpar_defaults()
{
  std::vector<double> dp(3);
  int n = dp.size();
  mgt_get_rpar_defaults(&dp[0], &n);
  return dp;
}
