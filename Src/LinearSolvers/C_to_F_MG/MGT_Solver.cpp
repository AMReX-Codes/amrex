#include <ParmParse.H>
#include <MGT_Solver.H>

bool  MGT_Solver::initialized;
int   MGT_Solver::def_nu_1;
int   MGT_Solver::def_nu_2;
int   MGT_Solver::def_nu_b;
int   MGT_Solver::def_nu_f;
int   MGT_Solver::def_gamma;
Real  MGT_Solver::def_omega;
int   MGT_Solver::def_maxiter;
int   MGT_Solver::def_maxiter_b;
int   MGT_Solver::def_verbose;
int   MGT_Solver::def_cg_verbose;
int   MGT_Solver::def_max_nlevel;
int   MGT_Solver::def_min_width;
int   MGT_Solver::def_smoother;
int   MGT_Solver::def_cycle;
int   MGT_Solver::def_usecg;
int   MGT_Solver::def_cg_solver;
int   MGT_Solver::def_bottom_solver;
Real  MGT_Solver::def_bottom_solver_eps;

typedef void (*mgt_get)(const int* lev, const int* n, double* uu, 
			const int* plo, const int* phi, 
			const int* lo, const int* hi);
typedef void (*mgt_set)(const int* lev, const int* n, const double* uu, 
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
		       bool nodal)
  :
  m_bd(bd), m_dmap(dmap), m_grids(grids), m_nodal(nodal)
{
   if (!initialized)
        initialize();

  BL_ASSERT( m_grids.size() == dmap.size() );
  m_nlevel = grids.size();
  int dm = BL_SPACEDIM;
  int i_nodal = (m_nodal)?1:0;

  mgt_alloc(&dm, &m_nlevel, &i_nodal);

  mgt_set_defaults(&def_nu_1,&def_nu_2,&def_nu_b,&def_nu_f,&def_gamma,&def_omega,
                   &def_maxiter,&def_maxiter_b,&def_bottom_solver,&def_bottom_solver_eps,
                   &def_verbose,&def_cg_verbose,&def_max_nlevel,&def_min_width,&def_cycle,&def_smoother);

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
      mgt_set_level(&lev, &nb, &dm, &lo[0], &hi[0], 
		    domain.loVect(), domain.hiVect(), &bc[0], pm, &pmap[0]);
    }

  mgt_finalize(dx);
}

void
MGT_Solver::initialize()
{
    initialized = true;

    mgt_init();

    mgt_get_defaults(&def_nu_1,&def_nu_2,&def_nu_b,&def_nu_f,&def_gamma,&def_omega,
                     &def_maxiter,&def_maxiter_b,&def_bottom_solver,
                     &def_verbose,&def_cg_verbose,&def_max_nlevel,&def_min_width,&def_cycle,&def_smoother);

    /* SET TO AGREE WITH MULTIGRID DEFAULT */
    def_min_width = 2;
    def_usecg = 1;
    def_cg_solver = 1;
    def_bottom_solver_eps = 0.0001;

    ParmParse pp("mg");

    pp.query("maxiter", def_maxiter);
    pp.query("maxiter_b", def_maxiter_b);
    pp.query("nu_1", def_nu_1);
    pp.query("nu_2", def_nu_2);
    pp.query("nu_b", def_nu_b);
    pp.query("nu_f", def_nu_f);
    pp.query("v"   , def_verbose);
    pp.query("usecg", def_usecg);
    pp.query("cg_solver", def_cg_solver);
    pp.query("rtol_b", def_bottom_solver_eps);
    pp.query("numLevelsMAX", def_max_nlevel);

/*
    pp.query("nu_0", def_nu_0);
    pp.query("bot_atol", def_atol_b);
    pp.query("smooth_on_cg_unstable", def_smooth_on_cg_unstable);
*/

    {
    ParmParse pp("cg");
    pp.query("v"   , def_cg_verbose);
    }

    if (def_usecg == 1) {
      if (def_cg_solver == 1) {
        def_bottom_solver = 1;
      } else if (def_cg_solver == 0) {
        def_bottom_solver = 2;
      }
    } else {
      def_bottom_solver = 0;
    }
}

void
MGT_Solver::set_coefficients(const MultiFab* aa[], const MultiFab* bb[][BL_SPACEDIM])
{
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      mgt_init_coeffs_lev(&lev);
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
	   mgt_set_cfa (&lev, &n, a->dataPtr(), alo, ahi, lo, hi);

	   const int* bxlo = b[0]->box().loVect();
	   const int* bxhi = b[0]->box().hiVect();
	   mgt_set_cfbx(&lev, &n, b[0]->dataPtr(), bxlo, bxhi, lo, hi);

	   const int* bylo = b[1]->box().loVect();
	   const int* byhi = b[1]->box().hiVect();
	   mgt_set_cfby(&lev, &n, b[1]->dataPtr(), bylo, byhi, lo, hi);

#if (BL_SPACEDIM == 3)
           const int* bzlo = b[2]->box().loVect();
  	   const int* bzhi = b[2]->box().hiVect();
  	   mgt_set_cfbz(&lev, &n, b[2]->dataPtr(), bzlo, bzhi, lo, hi);
#endif
	}
      mgt_finalize_stencil_lev(&lev, xa, xb, pxa, pxb);
    }
  mgt_finalize_stencil();
}

void 
MGT_Solver::solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol)
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
	  mgt_set_rh(&lev, &n, rd, rlo, rhi, lo, hi);

	  const FArrayBox& sol = (*(uu[lev]))[umfi];
	  const Real* sd = sol.dataPtr();
	  const int* slo = sol.box().loVect();
	  const int* shi = sol.box().hiVect();
	  mgt_set_uu(&lev, &n, sd, slo, shi, lo, hi);
	}
    }

  mgt_solve(tol,abs_tol);

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
	  mgt_get_uu(&lev, &n, sd, plo, phi, lo, hi);
	}
    }
}

MGT_Solver::~MGT_Solver()
{
  mgt_dealloc();
}
