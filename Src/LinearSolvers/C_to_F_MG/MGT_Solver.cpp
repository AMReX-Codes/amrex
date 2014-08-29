#include <ParmParse.H>
#include <MGT_Solver.H>
#include <ParallelDescriptor.H>

bool  MGT_Solver::initialized = false;
int   MGT_Solver::def_nu_1;
int   MGT_Solver::def_nu_2;
int   MGT_Solver::def_nu_b;
int   MGT_Solver::def_nu_f;
Real  MGT_Solver::def_max_L0_growth;
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
typedef void (*mgt_getni)(const int* lev, const int* n, double* uu, 
			const int* plo, const int* phi, 
			const int* lo, const int* hi, const int& nuu, const int&iuu);
typedef void (*mgt_get_ng)(const int* lev, const int* n, double* uu, 
			   const int* plo, const int* phi, 
   			   const int* lo, const int* hi, const int* ng);
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
			 const int* lo, const int* hi, const Real* r);
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
		          const int* lo, const int* hi, const Real* value);
#if BL_SPACEDIM == 1
mgt_get_ng  mgt_get_uu         = mgt_get_uu_1d;
mgt_set     mgt_set_uu         = mgt_set_uu_1d;
mgt_getni   mgt_get_pr         = mgt_get_pr_1d;
mgt_get     mgt_get_res        = mgt_get_res_1d;
mgt_get_dir mgt_get_gp         = mgt_get_gp_1d;
mgt_setni   mgt_set_pr         = mgt_set_pr_1d;
mgt_set     mgt_set_rh         = mgt_set_rh_1d;
mgt_set     mgt_set_cfa        = mgt_set_cfa_1d;
mgt_setn    mgt_set_cfa2       = mgt_set_cfa2_1d;
mgt_set_cf  mgt_set_cfbx       = mgt_set_cfbx_1d;
mgt_set_cfn mgt_set_cfbnx      = mgt_set_cfbnx_1d;
mgt_set_c   mgt_set_cfa_const  = mgt_set_cfa_1d_const;
mgt_set_c   mgt_set_cfbx_const = mgt_set_cfbx_1d_const;
mgt_set     mgt_set_cfs        = mgt_set_cfs_1d;
mgt_getni   mgt_get_vel        = mgt_get_vel_1d;
mgt_setni   mgt_set_vel        = mgt_set_vel_1d;
mgt_setr    mgt_add_rh_nodal   = mgt_add_rh_nodal_1d;
mgt_set     mgt_set_sync_msk   = mgt_set_sync_msk_1d;
mgt_set     mgt_set_vold       = mgt_set_vold_1d;
mgt_get     mgt_get_sync_res   = mgt_get_sync_res_1d;
mgt_set     mgt_set_rhcc_nodal = mgt_set_rhcc_nodal_1d;
#elif BL_SPACEDIM == 2
mgt_get_ng  mgt_get_uu         = mgt_get_uu_2d;
mgt_set     mgt_set_uu         = mgt_set_uu_2d;
mgt_getni   mgt_get_pr         = mgt_get_pr_2d;
mgt_get     mgt_get_res        = mgt_get_res_2d;
mgt_get_dir mgt_get_gp         = mgt_get_gp_2d;
mgt_setni   mgt_set_pr         = mgt_set_pr_2d;
mgt_set     mgt_set_rh         = mgt_set_rh_2d;
mgt_set     mgt_set_cfa        = mgt_set_cfa_2d;
mgt_setn    mgt_set_cfa2       = mgt_set_cfa2_2d;
mgt_set_c   mgt_set_cfa_const  = mgt_set_cfa_2d_const;
mgt_set_cf  mgt_set_cfbx       = mgt_set_cfbx_2d;
mgt_set_cfn mgt_set_cfbnx      = mgt_set_cfbnx_2d;
mgt_set_c   mgt_set_cfbx_const = mgt_set_cfbx_2d_const;
mgt_set_cf  mgt_set_cfby       = mgt_set_cfby_2d;
mgt_set_cfn mgt_set_cfbny      = mgt_set_cfbny_2d;
mgt_set_c   mgt_set_cfby_const = mgt_set_cfby_2d_const;
mgt_set     mgt_set_cfs        = mgt_set_cfs_2d;
mgt_getni   mgt_get_vel        = mgt_get_vel_2d;
mgt_setni   mgt_set_vel        = mgt_set_vel_2d;
mgt_setr    mgt_add_rh_nodal   = mgt_add_rh_nodal_2d;
mgt_set     mgt_set_sync_msk   = mgt_set_sync_msk_2d;
mgt_set     mgt_set_vold       = mgt_set_vold_2d;
mgt_get     mgt_get_sync_res   = mgt_get_sync_res_2d;
mgt_set     mgt_set_rhcc_nodal = mgt_set_rhcc_nodal_2d;
#elif BL_SPACEDIM == 3
mgt_get_ng  mgt_get_uu         = mgt_get_uu_3d;
mgt_set     mgt_set_uu         = mgt_set_uu_3d;
mgt_getni   mgt_get_pr         = mgt_get_pr_3d;
mgt_get     mgt_get_res        = mgt_get_res_3d;
mgt_get_dir mgt_get_gp         = mgt_get_gp_3d;
mgt_setni   mgt_set_pr         = mgt_set_pr_3d;
mgt_set     mgt_set_rh         = mgt_set_rh_3d;
mgt_set     mgt_set_cfa        = mgt_set_cfa_3d;
mgt_setn    mgt_set_cfa2       = mgt_set_cfa2_3d;
mgt_set_c   mgt_set_cfa_const  = mgt_set_cfa_3d_const;
mgt_set_cf  mgt_set_cfbx       = mgt_set_cfbx_3d;
mgt_set_cfn mgt_set_cfbnx      = mgt_set_cfbnx_3d;
mgt_set_c   mgt_set_cfbx_const = mgt_set_cfbx_3d_const;
mgt_set_cf  mgt_set_cfby       = mgt_set_cfby_3d;
mgt_set_cfn mgt_set_cfbny      = mgt_set_cfbny_3d;
mgt_set_c   mgt_set_cfby_const = mgt_set_cfby_3d_const;
mgt_set_cf  mgt_set_cfbz       = mgt_set_cfbz_3d;
mgt_set_cfn mgt_set_cfbnz      = mgt_set_cfbnz_3d;
mgt_set_c   mgt_set_cfbz_const = mgt_set_cfbz_3d_const;
mgt_set     mgt_set_cfs        = mgt_set_cfs_3d;
mgt_getni   mgt_get_vel        = mgt_get_vel_3d;
mgt_setni   mgt_set_vel        = mgt_set_vel_3d;
mgt_setr    mgt_add_rh_nodal   = mgt_add_rh_nodal_3d;
mgt_set     mgt_set_sync_msk   = mgt_set_sync_msk_3d;
mgt_set     mgt_set_vold       = mgt_set_vold_3d;
mgt_get     mgt_get_sync_res   = mgt_get_sync_res_3d;
mgt_set     mgt_set_rhcc_nodal = mgt_set_rhcc_nodal_3d;
#endif


//
// Constructing a solver for the following operator: 
//    (\alpha I - \beta \sum_i (1/b_i) \nabla a_i \nabla) \phi
//

MGT_Solver::MGT_Solver(const std::vector<Geometry>& geom, 
                       int* bc, 
		       const std::vector<BoxArray>& grids,
		       const std::vector<DistributionMapping>& dmap,
		       bool nodal,
		       int stencil_type,
		       bool _have_rhcc,
                       int nc,
                       int ncomp,
                       int verbose)
    :
    m_nlevel(grids.size()),
    m_grids(grids),
    m_nodal(nodal),
    have_rhcc(_have_rhcc)
{
    BL_ASSERT(geom.size()==m_nlevel);
    BL_ASSERT(dmap.size()==m_nlevel);
    Build(geom,bc,stencil_type,dmap,nc,ncomp,verbose);
}


void
MGT_Solver::Build(const std::vector<Geometry>& geom, 
                  int* bc, 
                  int stencil_type,
                  const std::vector<DistributionMapping>& dmap,
                  int nc,
                  int ncomp,
                  int verbose)
    
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
  int lverbose = (verbose > 0) ? verbose : def_verbose;

  if (m_nodal) {
    mgt_nodal_alloc(&dm, &m_nlevel, &stencil_type);
    mgt_set_nodal_defaults(&def_nu_1,&def_nu_2,&def_nu_b,&def_nu_f,
                           &def_maxiter,&def_maxiter_b,&def_bottom_solver,&def_bottom_solver_eps,
                           &lverbose,&def_cg_verbose,&def_max_nlevel,
                           &def_min_width,&def_cycle,&def_smoother,&stencil_type);
  } else {
    mgt_cc_alloc(&dm, &m_nlevel, &stencil_type);
    mgt_set_defaults(&def_nu_1,&def_nu_2,&def_nu_b,&def_nu_f,
                     &def_maxiter,&def_maxiter_b,&def_bottom_solver,&def_bottom_solver_eps,
                     &def_max_L0_growth,
                     &lverbose,&def_cg_verbose,&def_max_nlevel,
                     &def_min_width,&def_cycle,&def_smoother,&stencil_type);
  }

  Array<int> pm(dm);
  for ( int i = 0; i < dm; ++i ) 
    {
      pm[i] = geom[0].isPeriodic(i)? 1 : 0;
    }

  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      const Array<int>& pmap = dmap[lev].ProcessorMap();
      Box domain = geom[lev].Domain();

      int nb = m_grids[lev].size();
      std::vector<int> lo(nb*dm);
      std::vector<int> hi(nb*dm);

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

  std::vector<Real> dx(m_nlevel*dm);
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
    //mgt_finalize_n(&dx[0],&bc[0],&nc,&ncomp);
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
    BoxLib::ExecOnFinalize(MGT_Solver::Finalize);

    initialized = true;

    int comm = 0;

#ifdef BL_USE_MPI
    comm = MPI_Comm_c2f(ParallelDescriptor::Communicator());
#endif

    mgt_init(&comm);

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
    pp.query("cycle_type", def_cycle); // 1 -> F, 2 -> W, 3 -> V
    //
    // The C++ code usually sets CG solver type using cg.cg_solver.
    // We'll allow people to also use mg.cg_solver but pick up the former as well.
    //
    if (!pp.query("cg_solver", def_cg_solver))
    {
        ParmParse pp("cg");

        pp.query("cg_solver", def_cg_solver);
    }

/*
    pp.query("nu_0", def_nu_0);
    pp.query("bot_atol", def_atol_b);
    pp.query("smooth_on_cg_unstable", def_smooth_on_cg_unstable);
*/
    {
        ParmParse pp("cg");
        pp.query("v", def_cg_verbose);
    }

    {
        ParmParse pp("fabarray");
        int doit = 0;
        pp.query("do_alltoallv", doit);
        if (doit)
        {
            if (ParallelDescriptor::IOProcessor())
                std::cout << "Using Do_AllToAllV in fParallel code ...\n";
            mgt_use_alltoallv();
        }
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


void
MGT_Solver::set_mac_coefficients(const MultiFab* aa[], 
                                 const MultiFab* bb[][BL_SPACEDIM], 
                                 Array< Array<Real> >& xa,
                                 Array< Array<Real> >& xb)
{
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      mgt_init_coeffs_lev(&lev);
      double pxa[BL_SPACEDIM], pxb[BL_SPACEDIM];

      for ( int i = 0; i < BL_SPACEDIM; ++i ) 
	{
	  pxa[i] = pxb[i] = 0;
	}

      Real beta = 1.0;

//    NOTE: we only pass in aa here in order to get the validbox.
      for (MFIter amfi(*(aa[lev])); amfi.isValid(); ++amfi)
	{
	  const FArrayBox* b[BL_SPACEDIM];

	  int n = amfi.index();

	  for ( int i = 0; i < BL_SPACEDIM; ++i )
	    {
	      b[i] = &((*(bb[lev][i]))[amfi]);
	    }

 	   const int* lo = amfi.validbox().loVect();
	   const int* hi = amfi.validbox().hiVect();

	   const int* bxlo = b[0]->box().loVect();
	   const int* bxhi = b[0]->box().hiVect();
	   mgt_set_cfbx(&lev, &n, b[0]->dataPtr(), &beta, bxlo, bxhi, lo, hi);

#if (BL_SPACEDIM >= 2)
	   const int* bylo = b[1]->box().loVect();
	   const int* byhi = b[1]->box().hiVect();
	   mgt_set_cfby(&lev, &n, b[1]->dataPtr(), &beta, bylo, byhi, lo, hi);
#endif

#if (BL_SPACEDIM == 3)
           const int* bzlo = b[2]->box().loVect();
  	   const int* bzhi = b[2]->box().hiVect();
  	   mgt_set_cfbz(&lev, &n, b[2]->dataPtr(), &beta, bzlo, bzhi, lo, hi);
#endif
	}
      int dm = BL_SPACEDIM;
      mgt_finalize_stencil_lev(&lev, xa[lev].dataPtr(), xb[lev].dataPtr(), pxa, pxb, &dm);
    }
  mgt_finalize_stencil();
}

void
MGT_Solver::set_gravity_coefficients(Array< PArray<MultiFab> >& coeffs,
                                     Array< Array<Real> >& xa,
                                     Array< Array<Real> >& xb,
                                     int is_constant)
{

   if (is_constant == 1) 
   {
      set_const_gravity_coeffs(xa,xb);
   }
   else
   {

      double pxa[BL_SPACEDIM], pxb[BL_SPACEDIM];
      for ( int i = 0; i < BL_SPACEDIM; ++i ) 
         pxa[i] = pxb[i] = 0.;

      for ( int lev = 0; lev < m_nlevel; ++lev )
      {
      mgt_init_coeffs_lev(&lev);

//    NOTE: the sign convention is because the elliptic solver solves
//           (alpha MINUS del dot beta grad) phi = RHS
//           Here alpha is zero and we want to solve del dot grad phi = RHS,
//             which is equivalent to MINUS del dot (MINUS ONE) grad phi = RHS.
      Real value_zero =  0.0;
      Real value_one  = -1.0;

      for (MFIter mfi((coeffs[lev][0])); mfi.isValid(); ++mfi)
        {
           int n = mfi.index();
           const int* lo = m_grids[lev][n].loVect();
           const int* hi = m_grids[lev][n].hiVect();

           mgt_set_cfa_const (&lev, &n, lo, hi, &value_zero);
 
           const int* bxlo = coeffs[lev][0][n].box().loVect();
           const int* bxhi = coeffs[lev][0][n].box().hiVect();
           mgt_set_cfbx(&lev, &n, coeffs[lev][0][n].dataPtr(), &value_one, bxlo, bxhi, lo, hi);
 
#if (BL_SPACEDIM >= 2) 
           const int* bylo = coeffs[lev][1][n].box().loVect(); 
           const int* byhi = coeffs[lev][1][n].box().hiVect();
   	   mgt_set_cfby(&lev, &n, coeffs[lev][1][n].dataPtr(), &value_one, bylo, byhi, lo, hi);
#endif
 
#if (BL_SPACEDIM == 3)
           const int* bzlo = coeffs[lev][2][n].box().loVect();
           const int* bzhi = coeffs[lev][2][n].box().hiVect();
	   mgt_set_cfbz(&lev, &n, coeffs[lev][2][n].dataPtr(), &value_one, bzlo, bzhi, lo, hi);
#endif
        }

      int dm = BL_SPACEDIM;
      mgt_finalize_stencil_lev(&lev, xa[lev].dataPtr(), xb[lev].dataPtr(), pxa, pxb, &dm);
    }
    mgt_finalize_stencil();
   }
}

void
MGT_Solver::set_const_gravity_coeffs(Array< Array<Real> >& xa,
                                     Array< Array<Real> >& xb)
{
   double pxa[BL_SPACEDIM], pxb[BL_SPACEDIM];
   for ( int i = 0; i < BL_SPACEDIM; ++i ) 
      pxa[i] = pxb[i] = 0.;

   // NOTE: the sign convention is because the elliptic solver solves
   //        (alpha MINUS del dot beta grad) phi = RHS
   //        Here alpha is zero and we want to solve del dot grad phi = RHS,
   //        which is equivalent to MINUS del dot (MINUS ONE) grad phi = RHS.
   Real value_zero =  0.0;
   Real value_one  = -1.0;

   int dm = BL_SPACEDIM;
   for ( int lev = 0; lev < m_nlevel; ++lev )
   {
      mgt_finalize_const_stencil_lev(&lev, &value_zero, &value_one, 
                                     xa[lev].dataPtr(), xb[lev].dataPtr(), pxa, pxb, &dm);
   }
   mgt_finalize_stencil();
}

void
MGT_Solver::set_visc_coefficients(PArray<MultiFab>& aa, 
				  Array<PArray<MultiFab> >& bb, 
                                  const Real& beta, 
                                  Array< Array<Real> >& xa,
                                  Array< Array<Real> >& xb,
				  int index_order)
{
  for ( int lev = 0; lev < m_nlevel; ++lev )
  {
    mgt_init_coeffs_lev(&lev);

    double pxa[BL_SPACEDIM], pxb[BL_SPACEDIM];

    for ( int i = 0; i < BL_SPACEDIM; ++i ) 
      pxa[i] = pxb[i] = 0.;
  
    for (MFIter amfi(aa[lev]); amfi.isValid(); ++amfi)
    {      
      int n = amfi.index();

      const FArrayBox& a = aa[lev][amfi];
      const int* lo = amfi.validbox().loVect();
      const int* hi = amfi.validbox().hiVect();
      
      const int* alo = a.box().loVect();
      const int* ahi = a.box().hiVect();
      mgt_set_cfa (&lev, &n, a.dataPtr(), alo, ahi, lo, hi);

      const FArrayBox& bx = (index_order==0) ? bb[0][lev][amfi] : bb[lev][0][amfi];
      const int* bxlo = bx.box().loVect();
      const int* bxhi = bx.box().hiVect();
      mgt_set_cfbx(&lev, &n, bx.dataPtr(), &beta, bxlo, bxhi, lo, hi);

#if (BL_SPACEDIM >= 2)
      const FArrayBox& by = (index_order==0) ? bb[1][lev][amfi] : bb[lev][1][amfi];
      const int* bylo = by.box().loVect();
      const int* byhi = by.box().hiVect();
      mgt_set_cfby(&lev, &n, by.dataPtr(), &beta, bylo, byhi, lo, hi);
#endif

#if (BL_SPACEDIM == 3)
      const FArrayBox& bz = (index_order==0) ? bb[2][lev][amfi] : bb[lev][2][amfi];
      const int* bzlo = bz.box().loVect();
      const int* bzhi = bz.box().hiVect();
      mgt_set_cfbz(&lev, &n, bz.dataPtr(), &beta, bzlo, bzhi, lo, hi);
#endif
    }
    int dm = BL_SPACEDIM;
    mgt_finalize_stencil_lev(&lev, xa[lev].dataPtr(), xb[lev].dataPtr(), pxa, pxb, &dm);
  }
  mgt_finalize_stencil();
}

void
MGT_Solver::set_visc_coefficients(const MultiFab* aa[], const MultiFab* bb[][BL_SPACEDIM], 
                                  const Real& beta, 
                                  Array< Array<Real> >& xa,
                                  Array< Array<Real> >& xb)
{
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      mgt_init_coeffs_lev(&lev);

      double pxa[BL_SPACEDIM], pxb[BL_SPACEDIM];

      for ( int i = 0; i < BL_SPACEDIM; ++i ) 
	{
	  pxa[i] = pxb[i] = 0.;
	}

      for (MFIter amfi(*(aa[lev])); amfi.isValid(); ++amfi)
	{
	  const FArrayBox* a = &((*(aa[lev]))[amfi]);
	  const FArrayBox* b[BL_SPACEDIM];

	  int n = amfi.index();

	  for ( int i = 0; i < BL_SPACEDIM; ++i )
	    {
	      b[i] = &((*(bb[lev][i]))[amfi]);
	    }
 	   
 	   const int* lo = amfi.validbox().loVect();
	   const int* hi = amfi.validbox().hiVect();

	   const int* alo = a->box().loVect();
	   const int* ahi = a->box().hiVect();
	   mgt_set_cfa (&lev, &n, a->dataPtr(), alo, ahi, lo, hi);

	   const int* bxlo = b[0]->box().loVect();
	   const int* bxhi = b[0]->box().hiVect();
	   mgt_set_cfbx(&lev, &n, b[0]->dataPtr(), &beta, bxlo, bxhi, lo, hi);

#if (BL_SPACEDIM >= 2)
	   const int* bylo = b[1]->box().loVect();
	   const int* byhi = b[1]->box().hiVect();
	   mgt_set_cfby(&lev, &n, b[1]->dataPtr(), &beta, bylo, byhi, lo, hi);
#endif

#if (BL_SPACEDIM == 3)
           const int* bzlo = b[2]->box().loVect();
  	   const int* bzhi = b[2]->box().hiVect();
  	   mgt_set_cfbz(&lev, &n, b[2]->dataPtr(), &beta, bzlo, bzhi, lo, hi);
#endif
	}
      int dm = BL_SPACEDIM;
      mgt_finalize_stencil_lev(&lev, xa[lev].dataPtr(), xb[lev].dataPtr(), pxa, pxb, &dm);
    }
  mgt_finalize_stencil();
}

void
MGT_Solver::set_visc_coefficients(MultiFab* aa[], MultiFab* bb[][BL_SPACEDIM], 
                                  const Real& beta, 
                                  Array< Array<Real> >& xa,
                                  Array< Array<Real> >& xb)
{
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      mgt_init_coeffs_lev(&lev);

      double pxa[BL_SPACEDIM], pxb[BL_SPACEDIM];

      for ( int i = 0; i < BL_SPACEDIM; ++i ) 
	{
	  pxa[i] = pxb[i] = 0.;
	}

      for (MFIter amfi(*(aa[lev])); amfi.isValid(); ++amfi)
	{
	  const FArrayBox* a = &((*(aa[lev]))[amfi]);
	  const FArrayBox* b[BL_SPACEDIM];

	  int n = amfi.index();

	  for ( int i = 0; i < BL_SPACEDIM; ++i )
	    {
	      b[i] = &((*(bb[lev][i]))[amfi]);
	    }
 	   
 	   const int* lo = amfi.validbox().loVect();
	   const int* hi = amfi.validbox().hiVect();

	   const int* alo = a->box().loVect();
	   const int* ahi = a->box().hiVect();
	   mgt_set_cfa (&lev, &n, a->dataPtr(), alo, ahi, lo, hi);

	   const int* bxlo = b[0]->box().loVect();
	   const int* bxhi = b[0]->box().hiVect();
	   mgt_set_cfbx(&lev, &n, b[0]->dataPtr(), &beta, bxlo, bxhi, lo, hi);

#if (BL_SPACEDIM >= 2)
	   const int* bylo = b[1]->box().loVect();
	   const int* byhi = b[1]->box().hiVect();
	   mgt_set_cfby(&lev, &n, b[1]->dataPtr(), &beta, bylo, byhi, lo, hi);
#endif

#if (BL_SPACEDIM == 3)
           const int* bzlo = b[2]->box().loVect();
  	   const int* bzhi = b[2]->box().hiVect();
  	   mgt_set_cfbz(&lev, &n, b[2]->dataPtr(), &beta, bzlo, bzhi, lo, hi);
#endif
	}
      int dm = BL_SPACEDIM;
      mgt_finalize_stencil_lev(&lev, xa[lev].dataPtr(), xb[lev].dataPtr(), pxa, pxb, &dm);
    }
  mgt_finalize_stencil();
}


void MGT_Solver::set_maxorder(const int max_order)
{
  mgt_set_maxorder(&max_order);
}

void
MGT_Solver::set_porous_coefficients(PArray<MultiFab>& a1, 
				    PArray<MultiFab>& a2, 
                                    Array<PArray<MultiFab> >& bb, 
                                    const Real& beta, 
                                    Array< Array<Real> >& xa, 
				    Array< Array<Real> >& xb,
                                    int nc_opt)
{
  int nc = bb[0][0].nComp();
  for ( int lev = 0; lev < m_nlevel; ++lev )
  {
    mgt_init_mc_coeffs_lev(&lev,&nc,&nc_opt);
    double pxa[BL_SPACEDIM], pxb[BL_SPACEDIM];
    for ( int i = 0; i < BL_SPACEDIM; ++i ) 
      pxa[i] = pxb[i] = 0;
   
    for (MFIter amfi(a1[lev]); amfi.isValid(); ++amfi)
    {
      int n = amfi.index();

      const FArrayBox& af1 = a1[lev][amfi];
      const int* lo = amfi.validbox().loVect();
      const int* hi = amfi.validbox().hiVect();
	  
      const int* a1lo = af1.box().loVect();
      const int* a1hi = af1.box().hiVect();
      mgt_set_cfa (&lev, &n, af1.dataPtr(), a1lo, a1hi, lo, hi);

      if (nc_opt == 0)
      {
	FArrayBox& af2 = a2[lev][amfi];
	const int* a2lo = af2.box().loVect();
	const int* a2hi = af2.box().hiVect();
	mgt_set_cfa2 (&lev, &n, af2.dataPtr(), a2lo, a2hi, lo, hi, a2[lev].nComp());
      }
      const FArrayBox& bx = bb[0][lev][amfi];
      const int* bxlo = bx.box().loVect();
      const int* bxhi = bx.box().hiVect();
      mgt_set_cfbnx(&lev, &n, bx.dataPtr(), &beta, bxlo, bxhi, lo, hi, bx.nComp());

#if (BL_SPACEDIM >= 2)
      const FArrayBox& by = bb[1][lev][amfi];
      const int* bylo = by.box().loVect();
      const int* byhi = by.box().hiVect();
      mgt_set_cfbny(&lev, &n, by.dataPtr(), &beta, bylo, byhi, lo, hi, by.nComp());
#endif

#if (BL_SPACEDIM == 3)
      const FArrayBox& bz = bb[2][lev][amfi];
      const int* bzlo = bz.box().loVect();
      const int* bzhi =  bz.box().hiVect();
      mgt_set_cfbnz(&lev, &n, bz.dataPtr(), &beta, bzlo, bzhi, lo, hi, bz.nComp());
#endif
    }
    int dm = BL_SPACEDIM;
    mgt_mc_finalize_stencil_lev(&lev, xa[lev].dataPtr(), xb[lev].dataPtr(), 
				pxa, pxb, &dm, &nc_opt);
  }
  mgt_finalize_stencil();
}

void
MGT_Solver::set_porous_coefficients(const MultiFab* a1[], const MultiFab* a2[], 
                                    const MultiFab* bb[][BL_SPACEDIM], 
                                    const Real& beta, 
                                    Array< Array<Real> >& xa, 
				    Array< Array<Real> >& xb,
                                    int nc_opt)
{
  int nc = (*bb[0][0]).nComp();
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      mgt_init_mc_coeffs_lev(&lev,&nc,&nc_opt);
      double pxa[BL_SPACEDIM], pxb[BL_SPACEDIM];

      for ( int i = 0; i < BL_SPACEDIM; ++i ) 
	{
	  pxa[i] = pxb[i] = 0;
	}

      for (MFIter amfi(*(a1[lev])); amfi.isValid(); ++amfi)
	{
	  const FArrayBox* af1 = &((*(a1[lev]))[amfi]);
	  const FArrayBox* b[BL_SPACEDIM];

	  int n = amfi.index();

	  for ( int i = 0; i < BL_SPACEDIM; ++i )
	  {
	    b[i] = &((*(bb[lev][i]))[amfi]);
	  }

	  const int* lo = amfi.validbox().loVect();
	  const int* hi = amfi.validbox().hiVect();
	  
	  const int* a1lo = af1->box().loVect();
	  const int* a1hi = af1->box().hiVect();
	  mgt_set_cfa (&lev, &n, af1->dataPtr(), a1lo, a1hi, lo, hi);
	  
	  if (nc_opt == 0)
	  {
	    const FArrayBox* af2 = &((*(a2[lev]))[amfi]);
	    const int* a2lo = af2->box().loVect();
	    const int* a2hi = af2->box().hiVect();
	    mgt_set_cfa2 (&lev, &n, af2->dataPtr(), a2lo, a2hi, lo, hi, a2[0]->nComp());
	  }

	  const int* bxlo = b[0]->box().loVect();
	  const int* bxhi = b[0]->box().hiVect();
	  mgt_set_cfbnx(&lev, &n, b[0]->dataPtr(), &beta, bxlo, bxhi, lo, hi, b[0]->nComp());

#if (BL_SPACEDIM >= 2)
	  const int* bylo = b[1]->box().loVect();
	  const int* byhi = b[1]->box().hiVect();
	  mgt_set_cfbny(&lev, &n, b[1]->dataPtr(), &beta, bylo, byhi, lo, hi, b[1]->nComp());
#endif

#if (BL_SPACEDIM == 3)
	  const int* bzlo = b[2]->box().loVect();
	  const int* bzhi = b[2]->box().hiVect();
	  mgt_set_cfbnz(&lev, &n, b[2]->dataPtr(), &beta, bzlo, bzhi, lo, hi, b[2]->nComp());
#endif
	}
      int dm = BL_SPACEDIM;
      
      mgt_mc_finalize_stencil_lev(&lev, xa[lev].dataPtr(), xb[lev].dataPtr(), 
				  pxa, pxb, &dm, &nc_opt);
    }
  mgt_finalize_stencil();
}

void
MGT_Solver::set_porous_coefficients(MultiFab* a1[], const  MultiFab* a2[], 
                                    MultiFab* bb[][BL_SPACEDIM], 
                                    const Real& beta, 
				    Array< Array<Real> >& xa,
				    Array< Array<Real> >& xb,
				    int nc_opt)
{
  int nc = (*bb[0][0]).nComp();
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      mgt_init_mc_coeffs_lev(&lev,&nc,&nc_opt);
      double pxa[BL_SPACEDIM], pxb[BL_SPACEDIM];

      for ( int i = 0; i < BL_SPACEDIM; ++i ) 
	{
	  pxa[i] = pxb[i] = 0;
	}

      for (MFIter amfi(*(a1[lev])); amfi.isValid(); ++amfi)
	{
	  const FArrayBox* af1 = &((*(a1[lev]))[amfi]);
	  
	  const FArrayBox* b[BL_SPACEDIM];
	  for ( int i = 0; i < BL_SPACEDIM; ++i )
	    {
	      b[i] = &((*(bb[lev][i]))[amfi]);
	    }

 	   int n = amfi.index();
 	   const int* lo = amfi.validbox().loVect();
	   const int* hi = amfi.validbox().hiVect();

	   const int* a1lo = af1->box().loVect();
	   const int* a1hi = af1->box().hiVect();
	   mgt_set_cfa (&lev, &n, af1->dataPtr(), a1lo, a1hi, lo, hi);

	   if (nc_opt == 0)
	   {
	     const FArrayBox* af2 = &((*(a2[lev]))[amfi]);
	     const int* a2lo = af2->box().loVect();
	     const int* a2hi = af2->box().hiVect();
	     mgt_set_cfa2 (&lev, &n, af2->dataPtr(), a2lo, a2hi, lo, hi, a2[0]->nComp());
	   }

	   const int* bxlo = b[0]->box().loVect();
	   const int* bxhi = b[0]->box().hiVect();
	   mgt_set_cfbnx(&lev, &n, b[0]->dataPtr(), &beta, bxlo, bxhi, lo, hi, b[0]->nComp());


#if (BL_SPACEDIM >= 2)
	   const int* bylo = b[1]->box().loVect();
	   const int* byhi = b[1]->box().hiVect();
	   mgt_set_cfbny(&lev, &n, b[1]->dataPtr(), &beta, bylo, byhi, lo, hi, b[1]->nComp());
#endif

#if (BL_SPACEDIM == 3)
           const int* bzlo = b[2]->box().loVect();
  	   const int* bzhi = b[2]->box().hiVect();
  	   mgt_set_cfbnz(&lev, &n, b[2]->dataPtr(), &beta, bzlo, bzhi, lo, hi, b[2]->nComp());
#endif
	}
      int dm = BL_SPACEDIM;
      mgt_mc_finalize_stencil_lev(&lev, xa[lev].dataPtr(), xb[lev].dataPtr(), 
				  pxa, pxb, &dm, &nc_opt);
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

// *******************************************************************************
// These four do *not* take a status flag
// *******************************************************************************

void 
MGT_Solver::solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
                  const BndryData& bd, Real& final_resnorm)
{
    solve(uu,rh,tol,abs_tol,bd,0,final_resnorm);
}

void 
MGT_Solver::solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
                  const BndryData& bd, int need_grad_phi, Real& final_resnorm)
{
  BL_PROFILE("MGT_Solver::solve(1)");
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
  uu[lev]->FillBoundary();

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

  // Pass in the status flag from here so we can know whehter the 
  //      solver converged
  int status = 0;
  mgt_solve(tol,abs_tol,&need_grad_phi,&final_resnorm,&status);

  if (status != 0) 
     BoxLib::Error("Multigrid did not converge!");

  int ng = 0;
  if (need_grad_phi == 1) ng = 1;

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
	  mgt_get_uu(&lev, &n, sd, plo, phi, lo, hi, &ng);
	}
    }
}

void 
MGT_Solver::solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
                  int need_grad_phi, Real& final_resnorm)
{
  BL_PROFILE("MGT_Solver::solve(2)");
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

  // Pass in the status flag from here so we can know whehter the 
  //      solver converged
  int status = 0;
  mgt_solve(tol,abs_tol,&need_grad_phi,&final_resnorm,&status);

  if (status != 0) 
     BoxLib::Error("Multigrid did not converge!");

  int ng = 0;
  if (need_grad_phi == 1) ng = 1;

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
	  mgt_get_uu(&lev, &n, sd, plo, phi, lo, hi, &ng);
	}
    }
}

void 
MGT_Solver::solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
                  const BndryData bd[], int need_grad_phi, Real& final_resnorm)
{
  BL_PROFILE("MGT_Solver::solve(3)");
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      // Copy the boundary register values into the solution array to
      // be copied into F90
      
      for (OrientationIter oitr; oitr; ++oitr)
	{
	  const FabSet& fs = bd[lev].bndryValues(oitr());
	  for (MFIter umfi(*(uu[lev])); umfi.isValid(); ++umfi)
	    {
	      FArrayBox& dest = (*(uu[lev]))[umfi];
	      dest.copy(fs[umfi],fs[umfi].box());
	    }
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

  // Pass in the status flag from here so we can know whehter the 
  //      solver converged
  int status = 0;
  mgt_solve(tol,abs_tol,&need_grad_phi,&final_resnorm,&status);

  if (status != 0) 
     BoxLib::Error("Multigrid did not converge!");

  int ng = 0;
  if (need_grad_phi == 1) ng = 1;

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
	  mgt_get_uu(&lev, &n, sd, plo, phi, lo, hi, &ng);
	}
    }
}
 
// *******************************************************************************
// These two take a status flag
// *******************************************************************************

void 
MGT_Solver::solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
                  const BndryData& bd, int need_grad_phi, Real& final_resnorm, int& status)
{
  BL_PROFILE("MGT_Solver::solve(4)");
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
  uu[lev]->FillBoundary();

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
  mgt_solve(tol,abs_tol,&need_grad_phi,&final_resnorm,&status);
  if (status == 1) return;

  int ng = 0;
  if (need_grad_phi == 1) ng = 1;

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
	  mgt_get_uu(&lev, &n, sd, plo, phi, lo, hi, &ng);
	}
    }
}

void 
MGT_Solver::solve(MultiFab* uu[], MultiFab* rh[], const Real& tol, const Real& abs_tol,
                  int need_grad_phi, Real& final_resnorm,int& status)
{
  BL_PROFILE("MGT_Solver::solve(5)");
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

  mgt_solve(tol,abs_tol,&need_grad_phi,&final_resnorm,&status);
  if (status == 1) return;

  int ng = 0;
  if (need_grad_phi == 1) ng = 1;

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
	  mgt_get_uu(&lev, &n, sd, plo, phi, lo, hi, &ng);
	}
    }
}
 
// *******************************************************************************
// End of solve options
// *******************************************************************************

void 
MGT_Solver::applyop(MultiFab* uu[], MultiFab* res[], const BndryData& bd)
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
  uu[lev]->FillBoundary();

  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      for (MFIter umfi(*(uu[lev])); umfi.isValid(); ++umfi)
	{
	  int n = umfi.index();

	  const int* lo = umfi.validbox().loVect();
	  const int* hi = umfi.validbox().hiVect();

	  const FArrayBox& sol = (*(uu[lev]))[umfi];
	  const Real* sd = sol.dataPtr();
	  const int* slo = sol.box().loVect();
	  const int* shi = sol.box().hiVect();
	  mgt_set_uu(&lev, &n, sd, slo, shi, lo, hi);
	}
    }

  mgt_applyop();

  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      for (MFIter mfi(*(res[lev])); mfi.isValid(); ++mfi)
	{
	  FArrayBox& resfab = (*(res[lev]))[mfi];
	  Real* rd = resfab.dataPtr();
	  int n = mfi.index();
	  const int* lo = mfi.validbox().loVect();
	  const int* hi = mfi.validbox().hiVect();
	  const int* rlo = resfab.box().loVect();
	  const int* rhi = resfab.box().hiVect();
	  mgt_get_res(&lev, &n, rd, rlo, rhi, lo, hi);
	}
    }
}

void 
MGT_Solver::applyop(MultiFab* uu[], MultiFab* res[], const BndryData bd[])
{
  // Copy the boundary register values into the solution array to be
  // copied into F90
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      // Copy the boundary register values into the solution array to be copied into F90
      for (OrientationIter oitr; oitr; ++oitr)
	{
	  const FabSet& fs = bd[lev].bndryValues(oitr());
	  for (MFIter umfi(*(uu[lev])); umfi.isValid(); ++umfi)
	    {
	      FArrayBox& dest = (*(uu[lev]))[umfi];
	      dest.copy(fs[umfi],fs[umfi].box());
	    }
	}
      uu[lev]->FillBoundary();
    }

  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      for (MFIter umfi(*(uu[lev])); umfi.isValid(); ++umfi)
	{
	  int n = umfi.index();

	  const int* lo = umfi.validbox().loVect();
	  const int* hi = umfi.validbox().hiVect();

	  const FArrayBox& sol = (*(uu[lev]))[umfi];
	  const Real* sd = sol.dataPtr();
	  const int* slo = sol.box().loVect();
	  const int* shi = sol.box().hiVect();
	  mgt_set_uu(&lev, &n, sd, slo, shi, lo, hi);
	}
    }

  mgt_applyop();

  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      for (MFIter mfi(*(res[lev])); mfi.isValid(); ++mfi)
	{
	  FArrayBox& resfab = (*(res[lev]))[mfi];
	  Real* rd = resfab.dataPtr();
	  int n = mfi.index();
	  const int* lo = mfi.validbox().loVect();
	  const int* hi = mfi.validbox().hiVect();
	  const int* rlo = resfab.box().loVect();
	  const int* rhi = resfab.box().hiVect();
	  mgt_get_res(&lev, &n, rd, rlo, rhi, lo, hi);
	}
    }
}

void 
MGT_Solver::applybc(MultiFab* uu[], const BndryData& bd)
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
  uu[lev]->FillBoundary();

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
      mgt_set_uu(&lev, &n, sd, plo, phi, lo, hi);
    }
  }
}

void 
MGT_Solver::compute_residual(MultiFab* uu[], MultiFab* rh[], MultiFab* res[], const BndryData& bd)
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

  mgt_compute_residual();

  for ( int lev = 0; lev < m_nlevel; ++lev )
    {
      for (MFIter mfi(*(res[lev])); mfi.isValid(); ++mfi)
	{
	  FArrayBox& resfab = (*(res[lev]))[mfi];
	  Real* rd = resfab.dataPtr();
	  int n = mfi.index();
	  const int* lo = mfi.validbox().loVect();
	  const int* hi = mfi.validbox().hiVect();
	  const int* rlo = resfab.box().loVect();
	  const int* rhi = resfab.box().hiVect();
	  mgt_get_res(&lev, &n, rd, rlo, rhi, lo, hi);
	}
    }
}

void 
MGT_Solver::get_fluxes(int lev, PArray<MultiFab>& flux, const Real* dx)
{
  mgt_compute_flux(lev);

  for ( int dir = 0; dir < BL_SPACEDIM; ++dir )
    {
      for (MFIter mfi(flux[dir]); mfi.isValid(); ++mfi)
        {
          FArrayBox& gp = flux[dir][mfi];
          Real* gpd = gp.dataPtr();

          int n = mfi.index();
          const int* lo = mfi.validbox().loVect();
          const int* hi = mfi.validbox().hiVect();
          const int* gplo = gp.box().loVect();
          const int* gphi = gp.box().hiVect();

          mgt_get_gp(&lev, &dir, &n, gpd, gplo, gphi, lo, hi);
          gp.mult(dx[dir]);
        }
    }

  mgt_delete_flux(lev);
}

void 
MGT_Solver::get_fluxes(int lev, 
		       MultiFab* const* flux, 
		       const Real* dx)
{
  mgt_compute_flux(lev);

  for ( int dir = 0; dir < BL_SPACEDIM; ++dir )
    {
      for (MFIter mfi(*flux[dir]); mfi.isValid(); ++mfi)
        {
          FArrayBox& gp = (*flux[dir])[mfi];
          Real* gpd = gp.dataPtr();

          int n = mfi.index();
          const int* lo = mfi.validbox().loVect();
          const int* hi = mfi.validbox().hiVect();
          const int* gplo = gp.box().loVect();
          const int* gphi = gp.box().hiVect();

          mgt_get_gp(&lev, &dir, &n, gpd, gplo, gphi, lo, hi);
          gp.mult(dx[dir]);
        }
    }

  mgt_delete_flux(lev);
}

void 
MGT_Solver::nodal_project(MultiFab* p[], MultiFab* vel[], MultiFab* rhcc[], const PArray<MultiFab>& rhnd,
			  const Real& tol, const Real& abs_tol,
			  int* lo_inflow, int* hi_inflow)
{
  for ( int lev = 0; lev < m_nlevel; ++lev )
    {

      BL_ASSERT( (*p[lev]).nGrow() == 1 );
      //      BL_ASSERT( (*vel[lev]).nGrow() == 1 );

      const int ncomp_p = (*p[lev]).nComp();
      const int ncomp_vel = (*vel[lev]).nComp();

      for (MFIter vmfi(*(vel[lev])); vmfi.isValid(); ++vmfi)
	{
	  int n = vmfi.index();

	  const int* lo = vmfi.validbox().loVect();
	  const int* hi = vmfi.validbox().hiVect();

	  const FArrayBox& velfab = (*(vel[lev]))[vmfi];
	  const Real* vd = velfab.dataPtr();
	  const int* vlo = velfab.box().loVect();
	  const int* vhi = velfab.box().hiVect();
	  const int ivel = 0;
	  mgt_set_vel(&lev, &n, vd, vlo, vhi, lo, hi, ncomp_vel, ivel);
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
	  const int ip = 0;
	  mgt_set_pr(&lev, &n, sd, slo, shi, lo, hi, ncomp_p, ip);
	}
    }

  mgt_divu(lo_inflow, hi_inflow);

  Real r;
  Real rhmax = 0.0;

  bool added_nodal_rhs = false;

  for ( int lev = 0; lev < rhnd.size(); ++lev ) {
    if (rhnd.defined(lev)) {
      for (MFIter rmfi(rhnd[lev]); rmfi.isValid(); ++rmfi) {
        int n = rmfi.index();
        
        const int* lo = rmfi.validbox().loVect();
        const int* hi = rmfi.validbox().hiVect();
        
        const FArrayBox& rhsfab = (rhnd[lev])[rmfi];
        const Real* rhsd = rhsfab.dataPtr();
        const int* rhslo = rhsfab.box().loVect();
        const int* rhshi = rhsfab.box().hiVect();
        mgt_add_rh_nodal(&lev, &n, rhsd, rhslo, rhshi, lo, hi, &r);
        rhmax = std::max(rhmax,r);
      }
      added_nodal_rhs = true;
    }
  }

  if (added_nodal_rhs) 
  {
      ParallelDescriptor::ReduceRealMax(rhmax,ParallelDescriptor::IOProcessorNumber());
      if (ParallelDescriptor::IOProcessor())
          std::cout << " F90: Source norm after adding nodal RHS is " << rhmax << std::endl;
  }
  
  if (have_rhcc) {
    BL_ASSERT(rhcc[0] != 0);
    for ( int lev = 0; lev < m_nlevel; ++lev ) {
      for (MFIter mfi(*(rhcc[lev])); mfi.isValid(); ++mfi) {
	int n = mfi.index();

	  const int* lo = mfi.validbox().loVect();
	  const int* hi = mfi.validbox().hiVect();

	  const FArrayBox& rhccfab = (*(rhcc[lev]))[mfi];
	  const Real* p = rhccfab.dataPtr();
	  const int* plo = rhccfab.box().loVect();
	  const int* phi = rhccfab.box().hiVect();

	  mgt_set_rhcc_nodal(&lev, &n, p, plo, phi, lo, hi);
      }
    }

    mgt_add_divucc();
  }

  mgt_nodal_solve(tol,abs_tol);

  mgt_newu();

  for ( int lev = 0; lev < m_nlevel; ++lev )
    {

      const int ncomp_p = (*p[lev]).nComp();
      const int ncomp_vel = (*vel[lev]).nComp();

      for (MFIter pmfi(*(p[lev])); pmfi.isValid(); ++pmfi)
	{
	  FArrayBox& sol = (*(p[lev]))[pmfi];
	  Real* sd = sol.dataPtr();
	  int n = pmfi.index();
	  const int* lo = pmfi.validbox().loVect();
	  const int* hi = pmfi.validbox().hiVect();
	  const int* plo = sol.box().loVect();
	  const int* phi = sol.box().hiVect();
	  const int ip = 0;
	  mgt_get_pr(&lev, &n, sd, plo, phi, lo, hi, ncomp_p, ip);
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
	  const int ivel = 0;
	  mgt_get_vel(&lev, &n, vd, vlo, vhi, lo, hi, ncomp_vel, ivel);
	}
    }
}

void MGT_Solver::fill_sync_resid(MultiFab* sync_resid, const MultiFab& msk,
				 const MultiFab& vold, int isCoarse)
{
  mgt_alloc_nodal_sync();

  for (MFIter mfi(msk); mfi.isValid(); ++mfi) {
    int n = mfi.index();
    const FArrayBox& mskfab = msk[n];
    const Real* md = mskfab.dataPtr();
    const int* lo = mfi.validbox().loVect();
    const int* hi = mfi.validbox().hiVect();
    const int* plo = mskfab.box().loVect();
    const int* phi = mskfab.box().hiVect();
    const int lev = 0;
    mgt_set_sync_msk(&lev, &n, md, plo, phi, lo, hi);
  }

  for (MFIter mfi(vold); mfi.isValid(); ++mfi) {
    int n = mfi.index();
    const FArrayBox& vfab = vold[n];
    const Real* vd = vfab.dataPtr();
    const int* lo = mfi.validbox().loVect();
    const int* hi = mfi.validbox().hiVect();
    const int* plo = vfab.box().loVect();
    const int* phi = vfab.box().hiVect();
    const int lev = 0;
    mgt_set_vold(&lev, &n, vd, plo, phi, lo, hi);
  }

  if (isCoarse) {
    mgt_compute_sync_resid_crse();
  }
  else {
    mgt_compute_sync_resid_fine();
  }

  for (MFIter mfi(*sync_resid); mfi.isValid(); ++mfi) {
    int n = mfi.index();
    FArrayBox& sfab = (*sync_resid)[n];
    Real* sd = sfab.dataPtr();
    const int* lo = mfi.validbox().loVect();
    const int* hi = mfi.validbox().hiVect();
    const int* plo = sfab.box().loVect();
    const int* phi = sfab.box().hiVect();
    const int lev = 0;
    mgt_get_sync_res(&lev, &n, sd, plo, phi, lo, hi);
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
