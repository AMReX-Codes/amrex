#include <AMReX_ParmParse.H>
#include <AMReX_MGT_Solver.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

typedef void (*mgt_get)(const int* lev, const int* n, double* uu, 
			const int* plo, const int* phi, 
			const int* lo, const int* hi);
typedef void (*mgt_getni)(const int* lev, const int* n, double* uu, 
			const int* plo, const int* phi, 
			const int* lo, const int* hi, const int& nuu, const int&iuu);
typedef void (*mgt_get_dir)(const int* lev, const int* dir, const int* n, 
			    double* uu,
			    const int* plo, const int* phi, 
   			    const int* lo, const int* hi);
typedef void (*mgt_set)(const int* lev, const int* n, const double* uu, 
			const int* plo, const int* phi, 
			const int* lo, const int* hi);
typedef void (*mgt_setn)(const int* lev, const int* n, const double* uu, 
		 	 const int* plo, const int* phi, 
			 const int* lo, const int* hi, const int& nc);
typedef void (*mgt_setr)(const int* lev, const int* n, const double* uu, 
			 const int* plo, const int* phi, 
			 const int* lo, const int* hi, const amrex::Real* r);
typedef void (*mgt_setni)(const int* lev, const int* n, const double* uu, 
			  const int* plo, const int* phi, 
			  const int* lo, const int* hi, const int& nuu, const int& iuu);
typedef void (*mgt_set_cf)(const int* lev, const int* n, const double* uu, 
                           const double* b, 
                           const int* plo, const int* phi, 
                           const int* lo, const int* hi);
typedef void (*mgt_set_cfn)(const int* lev, const int* n, const double* uu, 
                            const double* b, 
                            const int* plo, const int* phi, 
                            const int* lo, const int* hi, const int& nc);
typedef void (*mgt_set_c)(const int* lev, const int* n, 
		          const int* lo, const int* hi, const amrex::Real* value);
#if BL_SPACEDIM == 1
mgt_get     mgt_get_uu         = mgt_get_uu_1d;
mgt_set     mgt_set_uu         = mgt_set_uu_1d;
mgt_getni   mgt_get_pr         = mgt_get_pr_1d;
mgt_get     mgt_get_res        = mgt_get_res_1d;
mgt_get_dir mgt_get_gp         = mgt_get_gp_1d;
mgt_setni   mgt_set_pr         = mgt_set_pr_1d;
mgt_set     mgt_set_rh         = mgt_set_rh_1d;
mgt_set     mgt_set_cfa        = mgt_set_cfa_1d;
mgt_setr    mgt_set_cfaa       = mgt_set_cfaa_1d;
mgt_set_cf  mgt_set_cfbx       = mgt_set_cfbx_1d;
mgt_set_cfn mgt_set_cfbnx      = mgt_set_cfbnx_1d;
mgt_set_c   mgt_set_cfa_const  = mgt_set_cfa_1d_const;
mgt_set     mgt_set_cfs        = mgt_set_cfs_1d;
mgt_getni   mgt_get_vel        = mgt_get_vel_1d;
mgt_setni   mgt_set_vel        = mgt_set_vel_1d;
mgt_setr    mgt_add_rh_nodal   = mgt_add_rh_nodal_1d;
mgt_set     mgt_set_rh_nodal   = mgt_set_rh_nodal_1d;
mgt_set     mgt_set_sync_msk   = mgt_set_sync_msk_1d;
mgt_set     mgt_set_vold       = mgt_set_vold_1d;
mgt_get     mgt_get_sync_res   = mgt_get_sync_res_1d;
mgt_set     mgt_set_rhcc_nodal = mgt_set_rhcc_nodal_1d;
#elif BL_SPACEDIM == 2
mgt_get     mgt_get_uu         = mgt_get_uu_2d;
mgt_set     mgt_set_uu         = mgt_set_uu_2d;
mgt_getni   mgt_get_pr         = mgt_get_pr_2d;
mgt_get     mgt_get_res        = mgt_get_res_2d;
mgt_get_dir mgt_get_gp         = mgt_get_gp_2d;
mgt_setni   mgt_set_pr         = mgt_set_pr_2d;
mgt_set     mgt_set_rh         = mgt_set_rh_2d;
mgt_set     mgt_set_cfa        = mgt_set_cfa_2d;
mgt_setr     mgt_set_cfaa      = mgt_set_cfaa_2d;
mgt_set_c   mgt_set_cfa_const  = mgt_set_cfa_2d_const;
mgt_set_cf  mgt_set_cfbx       = mgt_set_cfbx_2d;
mgt_set_cfn mgt_set_cfbnx      = mgt_set_cfbnx_2d;
mgt_set_cf  mgt_set_cfby       = mgt_set_cfby_2d;
mgt_set_cfn mgt_set_cfbny      = mgt_set_cfbny_2d;
mgt_set     mgt_set_cfs        = mgt_set_cfs_2d;
mgt_getni   mgt_get_vel        = mgt_get_vel_2d;
mgt_setni   mgt_set_vel        = mgt_set_vel_2d;
mgt_setr    mgt_add_rh_nodal   = mgt_add_rh_nodal_2d;
mgt_set     mgt_set_rh_nodal   = mgt_set_rh_nodal_2d;
mgt_set     mgt_set_sync_msk   = mgt_set_sync_msk_2d;
mgt_set     mgt_set_vold       = mgt_set_vold_2d;
mgt_get     mgt_get_sync_res   = mgt_get_sync_res_2d;
mgt_set     mgt_set_rhcc_nodal = mgt_set_rhcc_nodal_2d;
#elif BL_SPACEDIM == 3
mgt_get     mgt_get_uu         = mgt_get_uu_3d;
mgt_set     mgt_set_uu         = mgt_set_uu_3d;
mgt_getni   mgt_get_pr         = mgt_get_pr_3d;
mgt_get     mgt_get_res        = mgt_get_res_3d;
mgt_get_dir mgt_get_gp         = mgt_get_gp_3d;
mgt_setni   mgt_set_pr         = mgt_set_pr_3d;
mgt_set     mgt_set_rh         = mgt_set_rh_3d;
mgt_set     mgt_set_cfa        = mgt_set_cfa_3d;
mgt_setr     mgt_set_cfaa      = mgt_set_cfaa_3d;
mgt_set_c   mgt_set_cfa_const  = mgt_set_cfa_3d_const;
mgt_set_cf  mgt_set_cfbx       = mgt_set_cfbx_3d;
mgt_set_cfn mgt_set_cfbnx      = mgt_set_cfbnx_3d;
mgt_set_cf  mgt_set_cfby       = mgt_set_cfby_3d;
mgt_set_cfn mgt_set_cfbny      = mgt_set_cfbny_3d;
mgt_set_cf  mgt_set_cfbz       = mgt_set_cfbz_3d;
mgt_set_cfn mgt_set_cfbnz      = mgt_set_cfbnz_3d;
mgt_set     mgt_set_cfs        = mgt_set_cfs_3d;
mgt_getni   mgt_get_vel        = mgt_get_vel_3d;
mgt_setni   mgt_set_vel        = mgt_set_vel_3d;
mgt_setr    mgt_add_rh_nodal   = mgt_add_rh_nodal_3d;
mgt_set     mgt_set_rh_nodal   = mgt_set_rh_nodal_3d;
mgt_set     mgt_set_sync_msk   = mgt_set_sync_msk_3d;
mgt_set     mgt_set_vold       = mgt_set_vold_3d;
mgt_get     mgt_get_sync_res   = mgt_get_sync_res_3d;
mgt_set     mgt_set_rhcc_nodal = mgt_set_rhcc_nodal_3d;
#endif

namespace amrex {

bool  MGT_Solver::initialized = false;

int   MGT_Solver::def_maxiter;
int   MGT_Solver::def_maxiter_b;
int   MGT_Solver::def_bottom_solver;

int   MGT_Solver::def_nu_1;
int   MGT_Solver::def_nu_2;
int   MGT_Solver::def_nu_b;
int   MGT_Solver::def_nu_f;

int   MGT_Solver::def_verbose;
int   MGT_Solver::def_cg_verbose;

int   MGT_Solver::def_min_width;
int   MGT_Solver::def_max_nlevel;

int   MGT_Solver::def_cycle;
int   MGT_Solver::def_smoother;

int   MGT_Solver::def_usecg;
int   MGT_Solver::def_cg_solver;

Real  MGT_Solver::def_bottom_solver_eps;
Real  MGT_Solver::def_max_L0_growth;

//
// Constructing a solver for the following operator: 
//    (\alpha I - \beta \sum_i (1/b_i) \nabla a_i \nabla) \phi
//

MGT_Solver::MGT_Solver(const Vector<Geometry>& geom, 
                       int* bc, 
		       const Vector<BoxArray>& grids,
		       const Vector<DistributionMapping>& dmap,
		       bool nodal,
		       int stencil_type,
		       bool _have_rhcc,
                       int nc,
                       int ncomp,
                       int _verbose)
    :
    verbose(_verbose),
    m_nlevel(grids.size()),
    m_grids(grids),
    m_dmap(dmap),
    m_nodal(nodal),
    have_rhcc(_have_rhcc)
{
    BL_ASSERT(geom.size()==m_nlevel);
    BL_ASSERT(dmap.size()==m_nlevel);
    Build(geom,bc,stencil_type,dmap,nc,ncomp);
}


void
MGT_Solver::Build(const Vector<Geometry>& geom, 
                  int* bc, 
                  int stencil_type,
                  const Vector<DistributionMapping>& dmap,
                  int nc,
                  int ncomp)
    
{
   if (!initialized)
        initialize(m_nodal);

  BL_ASSERT(m_grids.size()==m_nlevel);
  BL_ASSERT(   dmap.size()==m_nlevel);
  int dm = BL_SPACEDIM;
  //
  // The default for "verbose" is false.
  // If it's true we use it since the user had to have set it somehow.
  // Otherwise we use def_verbose which is set generically using mg.v.
  //
  verbose = (verbose > 0) ? verbose : def_verbose;

  if (m_nodal) {
    mgt_nodal_alloc(&dm, &m_nlevel, &stencil_type);
    mgt_set_nodal_defaults(&def_nu_1,&def_nu_2,&def_nu_b,&def_nu_f,
                           &def_maxiter,&def_maxiter_b,&def_bottom_solver,&def_bottom_solver_eps,
                           &verbose,&def_cg_verbose,&def_max_nlevel,
                           &def_min_width,&def_cycle,&def_smoother,&stencil_type);
  } else {
    mgt_cc_alloc(&dm, &m_nlevel, &stencil_type);
    mgt_set_defaults(&def_nu_1,&def_nu_2,&def_nu_b,&def_nu_f,
                     &def_maxiter,&def_maxiter_b,&def_bottom_solver,&def_bottom_solver_eps,
                     &def_max_L0_growth,
                     &verbose,&def_cg_verbose,&def_max_nlevel,
                     &def_min_width,&def_cycle,&def_smoother,&stencil_type);
  }

  Vector<int> pm(dm);
  for ( int i = 0; i < dm; ++i ) 
    {
      pm[i] = geom[0].isPeriodic(i)? 1 : 0;
    }

  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      const Vector<int>& pmap = dmap[lev].ProcessorMap();
      Box domain = geom[lev].Domain();

      int nb = m_grids[lev].size();
      Vector<int> lo(nb*dm);
      Vector<int> hi(nb*dm);

      for ( int i = 0; i < nb; ++i )
      {
        for ( int j = 0; j < dm; ++j )
	{
	  lo[i + j*nb] = m_grids[lev][i].smallEnd(j);
	  hi[i + j*nb] = m_grids[lev][i].bigEnd(j);
	}
      }

      if (m_nodal) {
        mgt_set_nodal_level(&lev, &nb, &dm, &lo[0], &hi[0], 
  		            domain.loVect(), domain.hiVect(), pm.dataPtr(), &pmap[0]);
      } else {
        mgt_set_level(&lev, &nb, &dm, &lo[0], &hi[0], 
  		      domain.loVect(), domain.hiVect(), pm.dataPtr(), &pmap[0]);
      }
    }

  Vector<Real> dx(m_nlevel*dm);
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      for ( int j = 0; j < dm; ++j )
	{
	  dx[lev + j*m_nlevel] = geom[lev].CellSize()[j];
	}
    }

  if (m_nodal) {
    mgt_nodal_finalize(&dx[0],&bc[0]);
    if (have_rhcc) {
      mgt_alloc_rhcc_nodal();
    }
  } else {
    mgt_finalize(&dx[0],&bc[0]);
  }
}

void
MGT_Solver::Finalize()
{
    initialized = false;

    mgt_flush_copyassoc_cache();
}

void
MGT_Solver::FlushFortranOutput()
{
    mgt_flush_output();
}

void
MGT_Solver::initialize(bool nodal)
{
    amrex::ExecOnFinalize(MGT_Solver::Finalize);

    initialized = true;

    mgt_init();

    if (nodal) {
      mgt_get_nodal_defaults(&def_nu_1,&def_nu_2,&def_nu_b,&def_nu_f,
                             &def_maxiter,&def_maxiter_b,&def_bottom_solver,
                             &def_verbose,&def_cg_verbose,&def_max_nlevel,&def_min_width,&def_cycle,&def_smoother);
      def_nu_b = 2;
    } else {
      mgt_get_defaults(&def_nu_1,&def_nu_2,&def_nu_b,&def_nu_f,
                       &def_maxiter,&def_maxiter_b,&def_bottom_solver,&def_max_L0_growth,
                       &def_verbose,&def_cg_verbose,&def_max_nlevel,&def_min_width,&def_cycle,&def_smoother);
    }

    /* SET TO AGREE WITH MULTIGRID DEFAULT */
    def_min_width = 2;
    def_usecg = 1;
    def_cg_solver = 1;
    def_bottom_solver_eps = 0.0001;
    def_nu_f = 2;
    def_maxiter = 200;
    def_maxiter_b = 200;

    ParmParse pp("mg");

    pp.query("maxiter", def_maxiter);
    pp.query("maxiter_b", def_maxiter_b);
    pp.query("nu_1", def_nu_1);
    pp.query("nu_2", def_nu_2);
    pp.query("nu_b", def_nu_b);
    pp.query("nu_f", def_nu_f);
    pp.query("v"   , def_verbose);
    pp.query("usecg", def_usecg);

    pp.query("rtol_b", def_bottom_solver_eps);
    pp.query("numLevelsMAX", def_max_nlevel);
    pp.query("smoother", def_smoother);
    pp.query("cycle_type", def_cycle); // 1 -> F, 2 -> W, 3 -> V, 4 -> F+V
    //
    // The C++ code usually sets CG solver type using cg.cg_solver.
    // We'll allow people to also use mg.cg_solver but pick up the former as well.
    //
    if (!pp.query("cg_solver", def_cg_solver))
    {
        ParmParse ppcg("cg");

        ppcg.query("cg_solver", def_cg_solver);
    }

/*
    pp.query("nu_0", def_nu_0);
    pp.query("bot_atol", def_atol_b);
    pp.query("smooth_on_cg_unstable", def_smooth_on_cg_unstable);
*/
    {
        ParmParse ppcg("cg");
        ppcg.query("v", def_cg_verbose);
    }

    if (def_usecg == 1)
    {
        //
        // Translate from C++ -> F90 solver flag values.
        //
        if (def_cg_solver == 1)
        {
            //
            // BiCG
            //
            def_bottom_solver = 1;
        }
        else if (def_cg_solver == 0)
        {
            //
            // CG
            //
            def_bottom_solver = 2;
        }
        else if (def_cg_solver == 2)
        {
            //
            // CABiCG
            //
            def_bottom_solver = 3;
        }
    } else
    {
        //
        // Default to CABiCG.
        //
        def_bottom_solver = 3;
    }

    pp.query("bottom_solver", def_bottom_solver);
}


//
// (alpha * aa - beta * (del dot bb grad)) phi = RHS
// Here, aa is const one.
//
void
MGT_Solver::set_abeclap_coeffs (Real alpha,
				Real beta,
				const Vector< Vector<MultiFab*> >& bb,
				const Vector< Vector<Real> >& xa,
				const Vector< Vector<Real> >& xb)
{
    Vector<Real> pxa(BL_SPACEDIM, 0.0);
    Vector<Real> pxb(BL_SPACEDIM, 0.0);

    for ( int lev = 0; lev < m_nlevel; ++lev )
    {
	mgt_init_coeffs_lev(&lev);
    }

    for ( int lev = 0; lev < m_nlevel; ++lev )
    {
	set_cfa_const (alpha, lev);
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( int lev = 0; lev < m_nlevel; ++lev )
    {	
	for (int d=0; d<BL_SPACEDIM; ++d) {
	    set_cfb(*bb[lev][d], beta, lev, d);
	}
    }
	
    int dm = BL_SPACEDIM;
    for ( int lev = 0; lev < m_nlevel; ++lev )
    {
	mgt_finalize_stencil_lev(&lev, xa[lev].dataPtr(), xb[lev].dataPtr(), 
				 pxa.dataPtr(), pxb.dataPtr(), &dm);
    }

    mgt_finalize_stencil();
}

//
// (alpha * aa - beta * (del dot bb grad)) phi = RHS
// Here, alpha is one.
//
void
MGT_Solver::set_abeclap_coeffs (const Vector<MultiFab*>& aa,
				Real beta,
				const Vector< Vector<MultiFab*> >& bb,
				const Vector< Vector<Real> >& xa,
				const Vector< Vector<Real> >& xb)
{
    Vector<Real> pxa(BL_SPACEDIM, 0.0);
    Vector<Real> pxb(BL_SPACEDIM, 0.0);

    for ( int lev = 0; lev < m_nlevel; ++lev )
    {
	mgt_init_coeffs_lev(&lev);
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( int lev = 0; lev < m_nlevel; ++lev )
    {
	set_cfa(*aa[lev], lev);
	
	for (int d=0; d<BL_SPACEDIM; ++d) {
	    set_cfb(*bb[lev][d], beta, lev, d);
	}
    }
	
    int dm = BL_SPACEDIM;
    for ( int lev = 0; lev < m_nlevel; ++lev )
    {
	mgt_finalize_stencil_lev(&lev, xa[lev].dataPtr(), xb[lev].dataPtr(), 
				 pxa.dataPtr(), pxb.dataPtr(), &dm);
    }

    mgt_finalize_stencil();
}

//
// (alpha * aa - beta * (del dot bb grad)) phi = RHS
//
void
MGT_Solver::set_abeclap_coeffs (Real alpha,
				const Vector<MultiFab*>& aa,
				Real beta,
				const Vector< Vector<MultiFab*> >& bb,
				const Vector< Vector<Real> >& xa,
				const Vector< Vector<Real> >& xb)
{
    Vector<Real> pxa(BL_SPACEDIM, 0.0);
    Vector<Real> pxb(BL_SPACEDIM, 0.0);

    for ( int lev = 0; lev < m_nlevel; ++lev )
    {
	mgt_init_coeffs_lev(&lev);
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( int lev = 0; lev < m_nlevel; ++lev )
    {
	set_cfaa(*aa[lev], alpha, lev);
	
	for (int d=0; d<BL_SPACEDIM; ++d) {
	    set_cfb(*bb[lev][d], beta, lev, d);
	}
    }
	
    int dm = BL_SPACEDIM;
    for ( int lev = 0; lev < m_nlevel; ++lev )
    {
	mgt_finalize_stencil_lev(&lev, xa[lev].dataPtr(), xb[lev].dataPtr(), 
				 pxa.dataPtr(), pxb.dataPtr(), &dm);
    }

    mgt_finalize_stencil();
}

void
MGT_Solver::set_mac_coefficients(const Vector< Vector<MultiFab*> >& bb,
                                 const Vector< Vector<Real> >& xa,
                                 const Vector< Vector<Real> >& xb)
{
    Real alpha = 0.0;
    Real beta  = 1.0;
    set_abeclap_coeffs(alpha, beta, bb, xa, xb);
}

void
MGT_Solver::set_gravity_coefficients(const Vector< Vector<MultiFab*> >& bb,
                                     const Vector< Vector<Real> >& xa,
                                     const Vector< Vector<Real> >& xb)
{
    Real alpha =  0.0;
    Real beta  = -1.0;  // solving (del dot bb grad) phi = RHS
    set_abeclap_coeffs(alpha, beta, bb, xa, xb);
}

void
MGT_Solver::set_const_gravity_coeffs(const Vector< Vector<Real> >& xa,
                                     const Vector< Vector<Real> >& xb)
{
    Vector<Real> pxa(BL_SPACEDIM, 0.0);
    Vector<Real> pxb(BL_SPACEDIM, 0.0);

    Real alpha =  0.0;
    Real beta  = -1.0;  // solving (del dot grad) phi = RHS

    int dm = BL_SPACEDIM;
    for ( int lev = 0; lev < m_nlevel; ++lev )
    {
	mgt_finalize_const_stencil_lev(&lev, &alpha, &beta,
				       xa[lev].dataPtr(), xb[lev].dataPtr(), 
				       pxa.dataPtr(), pxb.dataPtr(), &dm);
    }

    mgt_finalize_stencil();
}

void
MGT_Solver::set_nodal_coefficients(const Vector<MultiFab*>& sig)
{
    for ( int lev = 0; lev < m_nlevel; ++lev ) {
	mgt_init_nodal_coeffs_lev(&lev);
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( int lev = 0; lev < m_nlevel; ++lev ) {
	for (MFIter mfi(*(sig[lev]),true); mfi.isValid(); ++mfi)
	{
	    const int n = mfi.LocalIndex();
	    const Box& bx = mfi.tilebox();
	    const FArrayBox& s = (*(sig[lev]))[mfi];
	    const Box& sbx = s.box();
	    mgt_set_cfs (&lev, &n, s.dataPtr(), sbx.loVect(), sbx.hiVect(), 
			 bx.loVect(), bx.hiVect());
	}
    }

    for ( int lev = 0; lev < m_nlevel; ++lev ) {
	mgt_finalize_nodal_stencil_lev(&lev);
    }

    mgt_finalize_nodal_stencil();
}

void
MGT_Solver::set_nodal_const_coefficients(Real val)
{
    for ( int lev = 0; lev < m_nlevel; ++lev ) {
	mgt_init_const_nodal_coeffs_lev(&lev,&val);
    }

    for ( int lev = 0; lev < m_nlevel; ++lev ) {
	mgt_finalize_nodal_stencil_lev(&lev);
    }

    mgt_finalize_nodal_stencil();
}

void MGT_Solver::set_maxorder(const int max_order)
{
  mgt_set_maxorder(&max_order);
}

void
MGT_Solver::solve(const Vector<MultiFab*>& uu, const Vector<MultiFab*>& rh, const BndryData& bd,
		  Real tol, Real abs_tol, int always_use_bnorm, 
		  Real& final_resnorm, int need_grad_phi)
{
    BL_PROFILE("MGT_Solver::solve()");

    // Copy the boundary register values into the solution array to be copied into F90
    int lev = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (OrientationIter oitr; oitr; ++oitr)
    {
	const FabSet& fs = bd.bndryValues(oitr());
	for (MFIter umfi(*(uu[lev])); umfi.isValid(); ++umfi)
	{
	    FArrayBox& dest = (*(uu[lev]))[umfi];
	    dest.copy(fs[umfi],fs[umfi].box());
	}
    }

#ifdef _OPENMP
#pragma omp parallel
#endif    
    for ( int ilev = 0; ilev < m_nlevel; ++ilev )
    {
	set_rh(*(rh[ilev]), ilev);
	set_uu(*(uu[ilev]), ilev);
    }
    
    // Pass in the status flag from here so we can know whether the 
    //      solver converged
    int status = 0;
    mgt_solve(tol,abs_tol,&need_grad_phi,&final_resnorm,&status,&always_use_bnorm);

    if (status != 0) 
	amrex::Error("Multigrid did not converge!");

    int ng = (need_grad_phi == 1) ? 1 : 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( int ilev = 0; ilev < m_nlevel; ++ilev )
    {
	get_uu(*(uu[ilev]), ilev, ng);
    }
}

void
MGT_Solver::solve_nodal(const Vector<MultiFab*>& uu, const Vector<MultiFab*>& rh,
                        Real tol, Real abs_tol)
{
    BL_PROFILE("MGT_Solver::solve_nodal()");

#ifdef _OPENMP
#pragma omp parallel
#endif    
    for ( int lev = 0; lev < m_nlevel; ++lev )
    {
	set_rh_nodal(*(rh[lev]), lev);
	set_uu_nodal(*(uu[lev]), lev);
    }
    
    mgt_nodal_solve(tol,abs_tol);
    
    int need_grad_phi = 0;
    int ng = (need_grad_phi == 1) ? 1 : 0;
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( int lev = 0; lev < m_nlevel; ++lev )
    {
	get_uu_nodal(*(uu[lev]), lev, ng);
    }
}

void 
MGT_Solver::applyop(const Vector<MultiFab*>& uu, const Vector<MultiFab*>& res, const BndryData& bd)
{
  // Copy the boundary register values into the solution array to be copied into F90
  int lev = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (OrientationIter oitr; oitr; ++oitr)
  {
      const FabSet& fs = bd.bndryValues(oitr());
      for (MFIter umfi(*(uu[lev])); umfi.isValid(); ++umfi)
      {
        FArrayBox& dest = (*(uu[lev]))[umfi];
        dest.copy(fs[umfi],fs[umfi].box());
      }
  }

#ifdef _OPENMP
#pragma omp parallel
#endif    
  for ( int ilev = 0; ilev < m_nlevel; ++ilev )
  {
      set_uu(*(uu[ilev]), ilev);
  }

  mgt_applyop();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( int ilev = 0; ilev < m_nlevel; ++ilev )
  {
      get_res(*(res[ilev]), ilev);
  }
}

void 
MGT_Solver::compute_residual(const Vector<MultiFab*>& uu, const Vector<MultiFab*>& rh,
			     const Vector<MultiFab*>& res, const BndryData& bd)
{
  // Copy the boundary register values into the solution array to be copied into F90
  int lev = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (OrientationIter oitr; oitr; ++oitr)
  {
      const FabSet& fs = bd.bndryValues(oitr());
      for (MFIter umfi(*(uu[lev])); umfi.isValid(); ++umfi)
      {
        FArrayBox& dest = (*(uu[lev]))[umfi];
        dest.copy(fs[umfi],fs[umfi].box());
      }
  }

#ifdef _OPENMP
#pragma omp parallel
#endif    
  for ( int ilev = 0; ilev < m_nlevel; ++ilev )
  {
      set_rh(*(rh[ilev]), ilev);
      set_uu(*(uu[ilev]), ilev);
  }

  mgt_compute_residual();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( int ilev = 0; ilev < m_nlevel; ++ilev )
  {
      get_res(*(res[ilev]), ilev);
  }
}

void 
MGT_Solver::get_fluxes(int lev, const Vector<MultiFab*>& flux, const Real* dx)
{
  mgt_compute_flux(lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( int dir = 0; dir < BL_SPACEDIM; ++dir )
  {
      get_gp(*flux[dir], lev, dir, dx[dir]);
  }

  mgt_delete_flux(lev);
}

void 
MGT_Solver::nodal_project(const Vector<MultiFab*>& p, const Vector<MultiFab*>& vel,
			  const Vector<MultiFab*>& rhcc, const Vector<MultiFab*>& rhnd,
			  const Real& tol, const Real& abs_tol,
			  int* lo_inflow, int* hi_inflow)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( int lev = 0; lev < m_nlevel; ++lev )
  {

      BL_ASSERT( (*p[lev]).nGrow() == 1 );
      //      BL_ASSERT( (*vel[lev]).nGrow() == 1 );

      const int ncomp_p = (*p[lev]).nComp();
      const int ncomp_vel = (*vel[lev]).nComp();

      for (MFIter mfi(*(vel[lev]), true); mfi.isValid(); ++mfi)
      {
	  const int n = mfi.LocalIndex();
	  const Box& ccbx = mfi.growntilebox(1);
	  const FArrayBox& v = (*(vel[lev]))[mfi];
	  const Box& vbox = v.box();
	  const int ivel = 0;
	  mgt_set_vel(&lev, &n, v.dataPtr(), vbox.loVect(), vbox.hiVect(), 
		      ccbx.loVect(), ccbx.hiVect(), ncomp_vel, ivel);

	  const Box& bx = mfi.grownnodaltilebox(-1,1);
	  const FArrayBox& pp = (*(p[lev]))[mfi];
	  const Box& pbox = pp.box();
	  const int ip = 0;
	  mgt_set_pr(&lev, &n, pp.dataPtr(), pbox.loVect(), pbox.hiVect(),
		     bx.loVect(), bx.hiVect(), ncomp_p, ip);
      }
  }

  mgt_divu(lo_inflow, hi_inflow);

  if (have_rhcc) BL_ASSERT(rhcc[0] != 0);

  Real rhmax = -1.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      Real rmax_this = -1.0;
      for ( int lev = 0; lev < rhnd.size(); ++lev )
      {
	  if (rhnd[lev]) {
	      for (MFIter mfi(*rhnd[lev], true); mfi.isValid(); ++mfi) 
	      {
		  const int n = mfi.LocalIndex();
		  const Box& bx = mfi.tilebox();
		  const FArrayBox& rfab = (*rhnd[lev])[mfi];
		  const Box& rbox = rfab.box();
		  mgt_add_rh_nodal(&lev, &n, rfab.dataPtr(), rbox.loVect(), rbox.hiVect(),
				   bx.loVect(), bx.hiVect(), &rmax_this);
	      }
	  }
      }

      if (have_rhcc) {
	  for ( int lev = 0; lev < m_nlevel; ++lev ) {
	      for (MFIter mfi(*(rhcc[lev]), true); mfi.isValid(); ++mfi) 
	      {
		  const int n = mfi.LocalIndex();
		  const Box& bx = mfi.tilebox();
		  const FArrayBox& rhccfab = (*(rhcc[lev]))[mfi];
		  const Box& rbox = rhccfab.box();
		  mgt_set_rhcc_nodal(&lev, &n, rhccfab.dataPtr(), rbox.loVect(), rbox.hiVect(),
				     bx.loVect(), bx.hiVect());
	      }
	  }
      }

      if (verbose > 0) {
#ifdef _OPENMP
#pragma omp critical (mgt_rhmax)
#endif
	  rhmax = std::max(rmax_this,rhmax);
      }
  }

  if (verbose > 0) {
      ParallelDescriptor::ReduceRealMax(rhmax,ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << " F90: Source norm after adding nodal RHS is " << rhmax << std::endl;
  }

  if (have_rhcc) {
    mgt_add_divucc();
  }

  mgt_nodal_solve(tol,abs_tol);

  mgt_newu();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( int lev = 0; lev < m_nlevel; ++lev )
  {

      const int ncomp_p = (*p[lev]).nComp();
      const int ncomp_vel = (*vel[lev]).nComp();

      for (MFIter mfi(*(vel[lev]), true); mfi.isValid(); ++mfi)
      {
	  const int n = mfi.LocalIndex();
	  const Box& ccbx = mfi.growntilebox(1);
	  FArrayBox& v = (*(vel[lev]))[mfi];
	  const Box& vbox = v.box();
	  const int ivel = 0;
	  mgt_get_vel(&lev, &n, v.dataPtr(), vbox.loVect(), vbox.hiVect(), 
		      ccbx.loVect(), ccbx.hiVect(), ncomp_vel, ivel);

	  const Box& bx = mfi.grownnodaltilebox(-1,1);
	  FArrayBox& pp = (*(p[lev]))[mfi];
	  const Box& pbox = pp.box();
	  const int ip = 0;
	  mgt_get_pr(&lev, &n, pp.dataPtr(), pbox.loVect(), pbox.hiVect(),
		     bx.loVect(), bx.hiVect(), ncomp_p, ip);
      }
  }
}

void MGT_Solver::fill_sync_resid(MultiFab& sync_resid, const MultiFab& msk,
				 const MultiFab& vold, int isCoarse)
{
  mgt_alloc_nodal_sync();

  const int lev = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(msk, true); mfi.isValid(); ++mfi) {
      const int n = mfi.LocalIndex();
      const Box& bx = mfi.growntilebox(1);

      const FArrayBox& mskfab = msk[mfi];
      const Box& mbox = mskfab.box();
      mgt_set_sync_msk(&lev, &n, mskfab.dataPtr(), mbox.loVect(), mbox.hiVect(),
		       bx.loVect(), bx.hiVect());

      const FArrayBox& vfab = vold[mfi];
      const Box& vbox = vfab.box();
      mgt_set_vold(&lev, &n, vfab.dataPtr(), vbox.loVect(), vbox.hiVect(),
		   bx.loVect(), bx.hiVect());
  }
  
  if (isCoarse) {
    mgt_compute_sync_resid_crse();
  }
  else {
    mgt_compute_sync_resid_fine();
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(sync_resid, true); mfi.isValid(); ++mfi) {
      const int n = mfi.LocalIndex();
      const Box& bx = mfi.tilebox();
      FArrayBox& sfab = sync_resid[mfi];
      const Box& sbx = sfab.box();
      mgt_get_sync_res(&lev, &n, sfab.dataPtr(), sbx.loVect(), sbx.hiVect(),
		       bx.loVect(), bx.hiVect());
  }

  mgt_dealloc_nodal_sync();
}

MGT_Solver::~MGT_Solver()
{
  if (m_nodal) {
    if (have_rhcc) {
      mgt_dealloc_rhcc_nodal();      
    }
    mgt_nodal_dealloc();
  } else {
    mgt_dealloc();
  }
}

void
MGT_Solver::set_cfa_const (Real alpha, int lev)
{
    if (alpha == 0.0) {
        // Nothing to do if alpha == 0, because cell coefficients are zero by default.
	return;
    } else {
	MFInfo info;
	info.SetAlloc(false);
	MultiFab cc(m_grids[lev], m_dmap[lev], 1, 0, info, FArrayBoxFactory()); // cell-centered MF      
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(cc, true); mfi.isValid(); ++mfi)
	{
	    int n = mfi.LocalIndex();
	    const Box& bx = mfi.tilebox();
	    mgt_set_cfa_const (&lev, &n, bx.loVect(), bx.hiVect(), &alpha);
	}
    }
}

void
MGT_Solver::set_cfa (const MultiFab& aa, int lev)
{
    // the caller has started OMP
    for (MFIter mfi(aa, true); mfi.isValid(); ++mfi)
    {
	const int n = mfi.LocalIndex();
	const Box& bx = mfi.tilebox();
	const FArrayBox& a = aa[mfi];
	const Box& abx = a.box();
	mgt_set_cfa (&lev, &n, a.dataPtr(),
		     abx.loVect(), abx.hiVect(),
		     bx.loVect(), bx.hiVect());
    }
}

void
MGT_Solver::set_cfaa (const MultiFab& aa, Real alpha, int lev)
{
    // the caller has started OMP
    for (MFIter mfi(aa, true); mfi.isValid(); ++mfi)
    {
	const int n = mfi.LocalIndex();
	const Box& bx = mfi.tilebox();
	const FArrayBox& a = aa[mfi];
	const Box& abx = a.box();
	mgt_set_cfaa (&lev, &n, a.dataPtr(),
		      abx.loVect(), abx.hiVect(),
		      bx.loVect(), bx.hiVect(),
		      &alpha);
    }
}

void 
MGT_Solver::set_cfb (const MultiFab& bb, Real beta, int lev, int dir)
{
    // the caller has started OMP
    for (MFIter mfi(bb,true); mfi.isValid(); ++mfi)
    {
	const int n = mfi.LocalIndex();
	const Box& bx = mfi.tilebox();
	const FArrayBox& fab = bb[mfi];
	const Box& fbox = fab.box();
	if (dir == 0) {
	    mgt_set_cfbx(&lev, &n, fab.dataPtr(), &beta, 
			 fbox.loVect(), fbox.hiVect(), 
			 bx.loVect(), bx.hiVect());
	}
#if (BL_SPACEDIM > 1)	
	else if (dir == 1) {
	    mgt_set_cfby(&lev, &n, fab.dataPtr(), &beta, 
			 fbox.loVect(), fbox.hiVect(), 
			 bx.loVect(), bx.hiVect());
	}
#if (BL_SPACEDIM == 3)	
	else {
	    mgt_set_cfbz(&lev, &n, fab.dataPtr(), &beta, 
			 fbox.loVect(), fbox.hiVect(), 
			 bx.loVect(), bx.hiVect());	    
	}
#endif
#endif
    }
}

void
MGT_Solver::set_rh (const MultiFab& mf, int lev)
{
    // the caller has started OMP
    for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
    {
	const int n = mfi.LocalIndex();
	const Box& bx = mfi.tilebox();
	const FArrayBox& fab = mf[mfi];
	const Box& fbx = fab.box();
	mgt_set_rh (&lev, &n, fab.dataPtr(),
		    fbx.loVect(), fbx.hiVect(),
		    bx.loVect(), bx.hiVect());
    }
}

void
MGT_Solver::set_uu (const MultiFab& mf, int lev)
{
    // the caller has started OMP
    for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
    {
	const int n = mfi.LocalIndex();
	const Box& bx = mfi.growntilebox(1);
	const FArrayBox& fab = mf[mfi];
	const Box& fbx = fab.box();
	mgt_set_uu (&lev, &n, fab.dataPtr(),
		    fbx.loVect(), fbx.hiVect(),
		    bx.loVect(), bx.hiVect());
    }
}

void
MGT_Solver::get_uu (MultiFab& mf, int lev, int ng)
{
    // the caller has started OMP
    for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
    {
	const int n = mfi.LocalIndex();
	const Box& bx = mfi.growntilebox(ng);
	FArrayBox& fab = mf[mfi];
	const Box& fbx = fab.box();
	mgt_get_uu (&lev, &n, fab.dataPtr(),
		    fbx.loVect(), fbx.hiVect(),
		    bx.loVect(), bx.hiVect());
    }
}

void
MGT_Solver::set_rh_nodal (const MultiFab& mf, int lev)
{
    // the caller has started OMP
    for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
    {
	const int n = mfi.LocalIndex();
	const Box& bx = mfi.tilebox();
	const FArrayBox& fab = mf[mfi];
	const Box& fbx = fab.box();
	mgt_set_rh_nodal (&lev, &n, fab.dataPtr(),
                          fbx.loVect(), fbx.hiVect(),
                          bx.loVect(), bx.hiVect());
    }
}

void
MGT_Solver::set_uu_nodal (const MultiFab& mf, int lev)
{
    // the caller has started OMP
    const int ncomp = mf.nComp();
    const int ip = 0;
    
    for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
    {
	const int n = mfi.LocalIndex();
	const Box& bx = mfi.growntilebox(1);
	const FArrayBox& fab = mf[mfi];
	const Box& fbx = fab.box();
	mgt_set_pr (&lev, &n, fab.dataPtr(),
                    fbx.loVect(), fbx.hiVect(),
                    bx.loVect(), bx.hiVect(), ncomp, ip);
    }
}

void
MGT_Solver::get_uu_nodal (MultiFab& mf, int lev, int ng)
{
    // the caller has started OMP
    const int ncomp = mf.nComp();
    const int ip = 0;

    for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
    {
	const int n = mfi.LocalIndex();
	const Box& bx = mfi.growntilebox(ng);
	FArrayBox& fab = mf[mfi];
	const Box& fbx = fab.box();
	mgt_get_pr (&lev, &n, fab.dataPtr(),
                    fbx.loVect(), fbx.hiVect(),
                    bx.loVect(), bx.hiVect(), ncomp, ip);
    }
}

void
MGT_Solver::get_res (MultiFab& mf, int lev)
{
    // the caller has started OMP
    for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
    {
	const int n = mfi.LocalIndex();
	const Box& bx = mfi.tilebox();
	FArrayBox& fab = mf[mfi];
	const Box& fbx = fab.box();
	mgt_get_res (&lev, &n, fab.dataPtr(),
		     fbx.loVect(), fbx.hiVect(),
		     bx.loVect(), bx.hiVect());
    }
}

void
MGT_Solver::get_gp (MultiFab& mf, int lev, int dir, Real dx)
{
    // the caller has started OMP
    for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
    {
	const int n = mfi.LocalIndex();
	const Box& bx = mfi.tilebox();
	FArrayBox& fab = mf[mfi];
	const Box& fbx = fab.box();
	mgt_get_gp (&lev, &dir, &n,
		    fab.dataPtr(), fbx.loVect(), fbx.hiVect(),
		    bx.loVect(), bx.hiVect());
	fab.mult(dx, bx);
    }
}

}
