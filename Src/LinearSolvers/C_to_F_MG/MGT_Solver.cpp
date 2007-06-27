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
int   MGT_Solver::stencil_type;

typedef void (*mgt_get)(const int* lev, const int* n, double* uu, 
			const int* plo, const int* phi, 
			const int* lo, const int* hi);
typedef void (*mgt_set)(const int* lev, const int* n, const double* uu, 
			const int* plo, const int* phi, 
			const int* lo, const int* hi);
typedef void (*mgt_set_cf)(const int* lev, const int* n, const double* uu, 
                           const double* b, 
                           const int* plo, const int* phi, 
                           const int* lo, const int* hi);
typedef void (*mgt_set_c)(const int* lev, const int* n, 
		          const int* lo, const int* hi, const Real* value);
#if BL_SPACEDIM == 2
mgt_get mgt_get_uu   = mgt_get_uu_2d;
mgt_set mgt_set_uu   = mgt_set_uu_2d;
mgt_get mgt_get_pr   = mgt_get_pr_2d;
mgt_set mgt_set_pr   = mgt_set_pr_2d;
mgt_set mgt_set_rh   = mgt_set_rh_2d;
mgt_set mgt_set_cfa  = mgt_set_cfa_2d;
mgt_set_cf mgt_set_cfbx = mgt_set_cfbx_2d;
mgt_set_cf mgt_set_cfby = mgt_set_cfby_2d;
mgt_set_c mgt_setc_cfa_const  = mgt_set_cfa_2d_const;
mgt_set_c mgt_setc_cfbx_const = mgt_set_cfbx_2d_const;
mgt_set_c mgt_setc_cfby_const = mgt_set_cfby_2d_const;
mgt_set mgt_set_cfs  = mgt_set_cfs_2d;
mgt_get mgt_get_vel  = mgt_get_vel_2d;
mgt_set mgt_set_vel  = mgt_set_vel_2d;
#elif BL_SPACEDIM == 3
mgt_get mgt_get_uu   = mgt_get_uu_3d;
mgt_set mgt_set_uu   = mgt_set_uu_3d;
mgt_get mgt_get_pr   = mgt_get_pr_3d;
mgt_set mgt_set_pr   = mgt_set_pr_3d;
mgt_set mgt_set_rh   = mgt_set_rh_3d;
mgt_set mgt_set_cfa  = mgt_set_cfa_3d;
mgt_set_cf mgt_set_cfbx = mgt_set_cfbx_3d;
mgt_set_cf mgt_set_cfby = mgt_set_cfby_3d;
mgt_set_cf mgt_set_cfbz = mgt_set_cfbz_3d;
mgt_set_c mgt_set_cfa_const  = mgt_set_cfa_3d_const;
mgt_set_c mgt_set_cfbx_const = mgt_set_cfbx_3d_const;
mgt_set_c mgt_set_cfby_const = mgt_set_cfby_3d_const;
mgt_set_c mgt_set_cfbz_const = mgt_set_cfbz_3d_const;
mgt_set_c mgt_set_cfb_constz = mgt_set_cfbz_3d_const;
mgt_set mgt_set_cfs  = mgt_set_cfs_3d;
mgt_get mgt_get_vel  = mgt_get_vel_2d;
mgt_set mgt_set_vel  = mgt_set_vel_3d;
#endif

MGT_Solver::MGT_Solver(const std::vector<Geometry>& geom, 
                       int* bc, 
		       const std::vector<BoxArray>& grids,
		       const std::vector<DistributionMapping>& dmap,
		       bool nodal)
  :
  m_dmap(dmap), m_grids(grids), m_nodal(nodal)
{
   if (!initialized)
        initialize(nodal);

  BL_ASSERT( m_grids.size() == dmap.size() );
  m_nlevel = grids.size();
  int dm = BL_SPACEDIM;
  int i_nodal = (m_nodal)?1:0;

  if (nodal) {
    mgt_nodal_alloc(&dm, &m_nlevel, &i_nodal, &stencil_type);
    mgt_set_nodal_defaults(&def_nu_1,&def_nu_2,&def_nu_b,&def_nu_f,&def_gamma,&def_omega,
                           &def_maxiter,&def_maxiter_b,&def_bottom_solver,&def_bottom_solver_eps,
                           &def_verbose,&def_cg_verbose,&def_max_nlevel,
                           &def_min_width,&def_cycle,&def_smoother,&stencil_type);
  } else {
    mgt_alloc(&dm, &m_nlevel, &i_nodal);
    mgt_set_defaults(&def_nu_1,&def_nu_2,&def_nu_b,&def_nu_f,&def_gamma,&def_omega,
                     &def_maxiter,&def_maxiter_b,&def_bottom_solver,&def_bottom_solver_eps,
                     &def_verbose,&def_cg_verbose,&def_max_nlevel,
                     &def_min_width,&def_cycle,&def_smoother,&stencil_type);
  }

  int nb = grids[0].size();
  std::vector<int> lo(nb*dm);
  std::vector<int> hi(nb*dm);
  int pm[dm];

  for ( int i = 0; i < dm; ++i ) 
    {
      pm[i] = geom[0].isPeriodic(i)? 1 : 0;
    }

  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      const Array<int>& pmap = dmap[lev].ProcessorMap();
      Box domain = geom[lev].Domain();

      for ( int i = 0; i < nb; ++i )
      {
        for ( int j = 0; j < dm; ++j )
	{
	  lo[i + j*nb] = m_grids[lev][i].smallEnd(j);
	  hi[i + j*nb] = m_grids[lev][i].bigEnd(j);
	}
      }

      if (nodal) {
        mgt_set_nodal_level(&lev, &nb, &dm, &lo[0], &hi[0], 
  		            domain.loVect(), domain.hiVect(), pm, &pmap[0]);
      } else {
        mgt_set_level(&lev, &nb, &dm, &lo[0], &hi[0], 
  		      domain.loVect(), domain.hiVect(), pm, &pmap[0]);
      }
    }

  std::vector<Real> dx(m_nlevel*dm);
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      for ( int j = 0; j < dm; ++j )
	{
	  dx[lev + j*m_nlevel] = geom[lev].CellSize()[j];
	}
    }

  if (nodal) {
    mgt_nodal_finalize(&dx[0],&bc[0]);
  } else {
    mgt_finalize(&dx[0],&bc[0]);
  }
}

void
MGT_Solver::initialize(bool nodal)
{
    initialized = true;

    mgt_init();

    if (nodal) {
      mgt_get_nodal_defaults(&def_nu_1,&def_nu_2,&def_nu_b,&def_nu_f,&def_gamma,&def_omega,
                             &def_maxiter,&def_maxiter_b,&def_bottom_solver,
                             &def_verbose,&def_cg_verbose,&def_max_nlevel,&def_min_width,&def_cycle,&def_smoother);
      def_nu_b = 2;
    } else {
      mgt_get_defaults(&def_nu_1,&def_nu_2,&def_nu_b,&def_nu_f,&def_gamma,&def_omega,
                       &def_maxiter,&def_maxiter_b,&def_bottom_solver,
                       &def_verbose,&def_cg_verbose,&def_max_nlevel,&def_min_width,&def_cycle,&def_smoother);
    }

    /* SET TO AGREE WITH MULTIGRID DEFAULT */
    def_min_width = 2;
    def_usecg = 1;
    def_cg_solver = 1;
    def_bottom_solver_eps = 0.0001;
    def_nu_f = 2;

    /* SET TO AGREE WITH MULTIGRID DEFAULT */
    stencil_type = 1;  /* 1: ST_CROSS */

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

    pp.query("stencil_type", stencil_type);

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
      def_bottom_solver = 3;
    }
}

void
MGT_Solver::set_mac_coefficients(const MultiFab* aa[], 
                                 const MultiFab* bb[][BL_SPACEDIM], const BndryData& bd)
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
          if (oitr().faceDir() == Orientation::low) {
            xa[dir] = bd.bndryLocs(oitr())[0];
          } else if (oitr().faceDir() == Orientation::high) {
            xb[dir] = bd.bndryLocs(oitr())[0];
          }
    
        }
 
      Real beta = 1.0;

//    NOTE: we only pass in aa here in order to get the validbox.
      for (MFIter amfi(*(aa[lev])); amfi.isValid(); ++amfi)
	{
	  const FArrayBox* b[BL_SPACEDIM];
	  for ( int i = 0; i < BL_SPACEDIM; ++i )
	    {
	      b[i] = &((*(bb[lev][i]))[amfi]);
	    }
 	   int n = amfi.index();
 	   const int* lo = amfi.validbox().loVect();
	   const int* hi = amfi.validbox().hiVect();

	   const int* bxlo = b[0]->box().loVect();
	   const int* bxhi = b[0]->box().hiVect();
	   mgt_set_cfbx(&lev, &n, b[0]->dataPtr(), &beta, bxlo, bxhi, lo, hi);

	   const int* bylo = b[1]->box().loVect();
	   const int* byhi = b[1]->box().hiVect();
	   mgt_set_cfby(&lev, &n, b[1]->dataPtr(), &beta, bylo, byhi, lo, hi);

#if (BL_SPACEDIM == 3)
           const int* bzlo = b[2]->box().loVect();
  	   const int* bzhi = b[2]->box().hiVect();
  	   mgt_set_cfbz(&lev, &n, b[2]->dataPtr(), &beta, bzlo, bzhi, lo, hi);
#endif
	}
      mgt_finalize_stencil_lev(&lev, xa, xb, pxa, pxb);
    }
  mgt_finalize_stencil();
}

void
MGT_Solver::set_gravity_coefficients(const BndryData& bd_crse)
{
  BoxArray grids(bd_crse.boxes());
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      mgt_init_coeffs_lev(&lev);
      double xa[BL_SPACEDIM], xb[BL_SPACEDIM];
      double pxa[BL_SPACEDIM], pxb[BL_SPACEDIM];

      for ( int i = 0; i < BL_SPACEDIM; ++i ) 
	{
	  pxa[i] = pxb[i] = 0;
	}

#if 0
      if (lev > 0) {
        for (int dir = 0; dir < BL_SPACEDIM; dir++) {
          xa[dir] = 0.5 * ref_ratio(lev-1,:)*mgt(lev)%dh(dir,mgt(lev)%nlevels)
          xb[dir] = 0.5 * ref_ratio(lev-1,:)*mgt(lev)%dh(dir,mgt(lev)%nlevels)
        }
      } else {
#endif
        for (OrientationIter oitr; oitr; ++oitr)
        {
          int dir  = oitr().coordDir();
          if (oitr().faceDir() == Orientation::low) {
            xa[dir] = bd_crse.bndryLocs(oitr())[0];
          } else if (oitr().faceDir() == Orientation::high) {
            xb[dir] = bd_crse.bndryLocs(oitr())[0];
          }
    
        }
#if 0
      }
#endif

      Real value_zero = 0.e0;
      Real value_one  = 1.e0;
      int ngrids = grids.size();
      for (int n = 0; n < ngrids; n++) 
        {
           const int* lo = grids[n].loVect();
           const int* hi = grids[n].hiVect();
           mgt_set_cfa_const (&lev, &n, lo, hi, &value_zero);
           mgt_set_cfbx_const(&lev, &n, lo, hi, &value_one);
           mgt_set_cfby_const(&lev, &n, lo, hi, &value_one);
#if (BL_SPACEDIM == 3)
           mgt_set_cfbz_const(&lev, &n, lo, hi, &value_one);
#endif
        }
      mgt_finalize_stencil_lev(&lev, xa, xb, pxa, pxb);
    }

 
  mgt_finalize_stencil();
}

void
MGT_Solver::set_visc_coefficients(const MultiFab* aa[], const MultiFab* bb[][BL_SPACEDIM], 
                                  const Real& beta, const BndryData& bd)
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
          if (oitr().faceDir() == Orientation::low) {
            xa[dir] = bd.bndryLocs(oitr())[0];
          } else if (oitr().faceDir() == Orientation::high) {
            xb[dir] = bd.bndryLocs(oitr())[0];
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
	   mgt_set_cfbx(&lev, &n, b[0]->dataPtr(), &beta, bxlo, bxhi, lo, hi);

	   const int* bylo = b[1]->box().loVect();
	   const int* byhi = b[1]->box().hiVect();
	   mgt_set_cfby(&lev, &n, b[1]->dataPtr(), &beta, bylo, byhi, lo, hi);

#if (BL_SPACEDIM == 3)
           const int* bzlo = b[2]->box().loVect();
  	   const int* bzhi = b[2]->box().hiVect();
  	   mgt_set_cfbz(&lev, &n, b[2]->dataPtr(), &beta, bzlo, bzhi, lo, hi);
#endif
	}
      mgt_finalize_stencil_lev(&lev, xa, xb, pxa, pxb);
    }
  mgt_finalize_stencil();
}

void
MGT_Solver::set_nodal_coefficients(const MultiFab* sig[])
{
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      mgt_init_nodal_coeffs_lev(&lev);

      for (MFIter smfi(*(sig[lev])); smfi.isValid(); ++smfi)
	{
	  const FArrayBox* s = &((*(sig[lev]))[smfi]);
 	  int n = smfi.index();
 	  const int* lo = smfi.validbox().loVect();
	  const int* hi = smfi.validbox().hiVect();

	  const int* slo = s->box().loVect();
	  const int* shi = s->box().hiVect();
	  mgt_set_cfs (&lev, &n, s->dataPtr(), slo, shi, lo, hi);
	}
      mgt_finalize_nodal_stencil_lev(&lev);
    }
  mgt_finalize_nodal_stencil();
}

void 
MGT_Solver::solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
                  const BndryData& bd)
{
  // Copy the boundary register values into the solution array to be copied into F90
  int lev = 0;
  for (OrientationIter oitr; oitr; ++oitr)
  {
      const FabSet& fs = bd.bndryValues(oitr());
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

void 
MGT_Solver::nodal_project(MultiFab* p[], MultiFab* vel[], const Real& tol, const Real& abs_tol)
{
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      for (MFIter vmfi(*(vel[lev])); vmfi.isValid(); ++vmfi)
	{
	  int n = vmfi.index();

	  const int* lo = vmfi.validbox().loVect();
	  const int* hi = vmfi.validbox().hiVect();

	  const FArrayBox& velfab = (*(vel[lev]))[vmfi];
	  const Real* vd = velfab.dataPtr();
	  const int* vlo = velfab.box().loVect();
	  const int* vhi = velfab.box().hiVect();
	  mgt_set_vel(&lev, &n, vd, vlo, vhi, lo, hi);
	}

      for (MFIter pmfi(*(p[lev])); pmfi.isValid(); ++pmfi)
	{
	  int n = pmfi.index();

	  const int* lo = pmfi.validbox().loVect();
	  const int* hi = pmfi.validbox().hiVect();

	  const FArrayBox& sol = (*(p[lev]))[pmfi];
	  const Real* sd = sol.dataPtr();
	  const int* slo = sol.box().loVect();
	  const int* shi = sol.box().hiVect();
	  mgt_set_pr(&lev, &n, sd, slo, shi, lo, hi);
	}
    }

  mgt_divu();

  mgt_nodal_solve(tol,abs_tol);

  mgt_newu();

  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      for (MFIter pmfi(*(p[lev])); pmfi.isValid(); ++pmfi)
	{
	  FArrayBox& sol = (*(p[lev]))[pmfi];
	  Real* sd = sol.dataPtr();
	  int n = pmfi.index();
	  const int* lo = pmfi.validbox().loVect();
	  const int* hi = pmfi.validbox().hiVect();
	  const int* plo = sol.box().loVect();
	  const int* phi = sol.box().hiVect();
	  mgt_get_pr(&lev, &n, sd, plo, phi, lo, hi);
	}

      for (MFIter vmfi(*(vel[lev])); vmfi.isValid(); ++vmfi)
	{
	  FArrayBox& velfab = (*(vel[lev]))[vmfi];
	  Real* vd = velfab.dataPtr();
	  int n = vmfi.index();
	  const int* lo = vmfi.validbox().loVect();
	  const int* hi = vmfi.validbox().hiVect();
	  const int* vlo = velfab.box().loVect();
	  const int* vhi = velfab.box().hiVect();
	  mgt_get_vel(&lev, &n, vd, vlo, vhi, lo, hi);
	}
    }
}

MGT_Solver::~MGT_Solver()
{
  if (m_nodal) {
    mgt_nodal_dealloc();
  } else {
    mgt_dealloc();
  }
}
