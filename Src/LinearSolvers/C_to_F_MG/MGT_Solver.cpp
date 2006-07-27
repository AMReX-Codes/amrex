#include <MGT_Solver.H>

typedef void (*mgt_get)(const int* mgt, const int* lev, const int* n, double* uu, 
			const int* plo, const int* phi, 
			const int* lo, const int* hi);
typedef void (*mgt_set)(const int* mgt, const int* lev, const int* n, const double* uu, 
			const int* plo, const int* phi, 
			const int* lo, const int* hi);
#if BL_SPACEDIM == 1
mgt_get mgt_get_uu   = mgt_get_uu_1d;
mgt_set mgt_set_uu   = mgt_set_uu_1d;
mgt_get mgt_get_rh   = mgt_get_rh_1d;
mgt_set mgt_set_rh   = mgt_set_rh_1d;
mgt_set mgt_set_cfa  = mgt_set_cfa_1d;
mgt_set mgt_set_cfbx = mgt_set_cfbx_1d;
#elif BL_SPACEDIM == 2
mgt_get mgt_get_uu   = mgt_get_uu_2d;
mgt_set mgt_set_uu   = mgt_set_uu_2d;
mgt_get mgt_get_rh   = mgt_get_rh_2d;
mgt_set mgt_set_rh   = mgt_set_rh_2d;
mgt_set mgt_set_cfa  = mgt_set_cfa_2d;
mgt_set mgt_set_cfbx = mgt_set_cfbx_2d;
mgt_set mgt_set_cfby = mgt_set_cfby_2d;
#elif BL_SPACEDIM == 3
mgt_get mgt_get_uu   = mgt_get_uu_3d;
mgt_set mgt_set_uu   = mgt_set_uu_3d;
mgt_get mgt_get_rh   = mgt_get_rh_3d;
mgt_set mgt_set_rh   = mgt_set_rh_3d;
mgt_set mgt_set_cfa  = mgt_set_cfa_3d;
mgt_set mgt_set_cfbx = mgt_set_cfbx_3d;
mgt_set mgt_set_cfby = mgt_set_cfby_3d;
mgt_set mgt_set_cfbz = mgt_set_cfbz_3d;
#endif
  
MGT_Solver::MGT_Solver(const BndryData& bd, const BCRec& phys_bc, const double* dx, 
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
	  // bc[i*2 + 0] = MGT_BC_DIR;	// FIXME: Hardware DIRICHET
	  // bc[i*2 + 1] = MGT_BC_DIR;
	  bc[i*2 + 0] = phys_bc.lo(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
	  bc[i*2 + 1] = phys_bc.hi(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
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
  mgt_finalize(&m_mgt, &n_ipar, ipar_p, &n_rpar, rpar_p, dx);
}

void
MGT_Solver::set_coefficients(const MultiFab* aa[], const MultiFab* bb[][BL_SPACEDIM])
{
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      mgt_init_coeffs_lev(&m_mgt, &lev);
      double xa[BL_SPACEDIM], xb[BL_SPACEDIM];
      double pxa[BL_SPACEDIM], pxb[BL_SPACEDIM];

      for ( int i = 0; i < BL_SPACEDIM; ++i ) 
	{
	  pxa[i] = pxb[i] = 0;
	}

      for (OrientationIter oitr; oitr; ++oitr)
        {
          int dir  = oitr().coordDir();
          int hilo = oitr().faceDir();
          if (oitr().faceDir() == Orientation::low) {
            xa[dir] = m_bd.bndryLocs(oitr())[0];
          } else if (oitr().faceDir() == Orientation::high) {
            xb[dir] = m_bd.bndryLocs(oitr())[0];
          }
    
        }

      for (MFIter amfi(*(aa[lev])); amfi.isValid(); ++amfi)
	{
	  const FArrayBox* a = &((*(aa[lev]))[amfi]);
	  const FArrayBox* b[BL_SPACEDIM];
	  for ( int i = 0; i < BL_SPACEDIM; ++i )
	    {
	      b[i] = &((*(bb[lev][i]))[amfi]);
	    }
 	   int n = amfi.index();
 	   const int* lo = amfi.validbox().loVect();
	   const int* hi = amfi.validbox().hiVect();

	   const int* alo = a->box().loVect();
	   const int* ahi = a->box().hiVect();
	   mgt_set_cfa (&m_mgt, &lev, &n, a->dataPtr(), alo, ahi, lo, hi);

	   const int* bxlo = b[0]->box().loVect();
	   const int* bxhi = b[0]->box().hiVect();
	   mgt_set_cfbx(&m_mgt, &lev, &n, b[0]->dataPtr(), bxlo, bxhi, lo, hi);

	   const int* bylo = b[1]->box().loVect();
	   const int* byhi = b[1]->box().hiVect();
	   mgt_set_cfby(&m_mgt, &lev, &n, b[1]->dataPtr(), bylo, byhi, lo, hi);
	}
      mgt_finalize_stencil_lev(&m_mgt, &lev, xa, xb, pxa, pxb);
    }
  mgt_finalize_stencil(&m_mgt);
}

void 
MGT_Solver::solve(MultiFab* uu[], MultiFab* rh[])
{
  // Copy the boundary register values into the solution array to be copied into F90
  int lev = 0;
  for (OrientationIter oitr; oitr; ++oitr)
  {
      const FabSet& fs = m_bd.bndryValues(oitr());
      for (MFIter umfi(*(uu[lev])); umfi.isValid(); ++umfi)
      {
        FArrayBox& dest = (*(uu[lev]))[umfi];
        dest.copy(fs[umfi],fs[umfi].box());
      }
  }

  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      for (MFIter umfi(*(uu[lev])); umfi.isValid(); ++umfi)
	{
	  int n = umfi.index();

	  const int* lo = umfi.validbox().loVect();
	  const int* hi = umfi.validbox().hiVect();

	  const FArrayBox& rhs = (*(rh[lev]))[umfi];
	  const Real* rd = rhs.dataPtr();
	  const int* rlo = rhs.box().loVect();
	  const int* rhi = rhs.box().hiVect();
	  mgt_set_rh(&m_mgt, &lev, &n, rd, rlo, rhi, lo, hi);

	  const FArrayBox& sol = (*(uu[lev]))[umfi];
	  const Real* sd = sol.dataPtr();
	  const int* slo = sol.box().loVect();
	  const int* shi = sol.box().hiVect();
	  mgt_set_uu(&m_mgt, &lev, &n, sd, slo, shi, lo, hi);
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
