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

void 
MGT_Solver::solve(MultiFab& uu, MultiFab& rh)
{
  int lev = 1;
  for (MFIter umfi(uu); umfi.isValid(); ++umfi)
    {
      const FArrayBox& rhs = rh[umfi];
      const FArrayBox& sol = uu[umfi];
      const Real* rd = rhs.dataPtr();
      const Real* sd = sol.dataPtr();
      int n = umfi.index() + 1; // Fortran Index Correction
      const int* lo = umfi.validbox().loVect();
      const int* hi = umfi.validbox().hiVect();
      const int* plo = rhs.box().loVect();
      const int* phi = rhs.box().hiVect();
      const int* splo = sol.box().loVect();
      const int* sphi = sol.box().hiVect();
      mgt_set_rh(&m_mgt, &lev, &n, rd, plo, phi, lo, hi);
      mgt_set_uu(&m_mgt, &lev, &n, sd, plo, phi, lo, hi);
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
  for (MFIter umfi(uu); umfi.isValid(); ++umfi)
    {
      FArrayBox& sol = uu[umfi];
      Real* sd = sol.dataPtr();
      int n = umfi.index() + 1; // Fortran Index Correction
      const int* lo = umfi.validbox().loVect();
      const int* hi = umfi.validbox().hiVect();
      const int* plo = sol.box().loVect();
      const int* phi = sol.box().hiVect();
      mgt_get_uu(&m_mgt, &lev, &n, sd, plo, phi, lo, hi);
    }
}

MGT_Solver::~MGT_Solver()
{
  mgt_dealloc(&m_mgt);
}
