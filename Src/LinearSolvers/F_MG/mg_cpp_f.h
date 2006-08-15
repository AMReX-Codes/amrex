#ifndef _MG_CPP_F_H_
#define _MG_CPP_F_H_

#if defined(BL_FORT_USE_UPPERCASE)
#define mgt_init                  MGT_INIT
#define mgt_alloc                 MGT_ALLOC
#define mgt_nodal_alloc           MGT_NODAL_ALLOC
#define mgt_set_level             MGT_SET_LEVEL
#define mgt_set_nodal_level       MGET_SET_NODAL_LEVEL
#define mgt_finalize              MGT_FINALIZE
#define mgt_nodal_finalize        MGT_NODAL_FINALIZE
#define mgt_init_coeffs_lev       MGT_INIT_COEFFS_LEV
#define mgt_init_nodal_coeffs_lev MGT_INIT_NODAL_COEFFS_LEV
#define mgt_finalize_stencil      MGT_FINALIZE_STENCIL
#define mgt_finalize_nodal_stencil MGT_FINALIZE_NODAL_STENCIL
#define mgt_finalize_stencil_lev  MGT_FINALIZE_STENCIL_LEV
#define mgt_finalize_nodal_stencil_lev  MGT_FINALIZE_NODAL_STENCIL_LEV
#define mgt_dealloc               MGT_DEALLOC
#define mgt_nodal_dealloc         MGT_NODAL_DEALLOC
#define mgt_solve                 MGT_SOLVE
#define mgt_nodal_solve           MGT_NODAL_SOLVE
#define mgt_divu                  MGT_DIVU
#define mgt_newu                  MGT_NEWU
#define mgt_set_defaults          MGT_SET_DEFAULTS
#define mgt_set_nodal_defaults    MGT_SET_NODAL_DEFAULTS
#define mgt_get_defaults          MGT_GET_DEFAULTS
#define mgt_get_nodal_defaults    MGT_GET_NODAL_DEFAULTS

#define mgt_set_rh_2d             MGT_SET_RH_2D
#define mgt_set_rh_3d             MGT_SET_RH_3D

#define mgt_set_uu_2d             MGT_SET_UU_2D
#define mgt_get_uu_2d             MGT_GET_UU_2D
#define mgt_set_uu_3d             MGT_SET_UU_3D
#define mgt_get_uu_3d             MGT_GET_UU_3D

#define mgt_set_pr_2d             MGT_SET_PR_2D
#define mgt_get_pr_2d             MGT_GET_PR_2D
#define mgt_set_pr_3d             MGT_SET_PR_3D
#define mgt_get_pr_3d             MGT_GET_PR_3D

#define mgt_set_cfa_2d            MGT_SET_CFA_2D
#define mgt_set_cfbx_2d           MGT_SET_CFBX_2D
#define mgt_set_cfby_2d           MGT_SET_CFBY_2D
#define mgt_set_cfa_3d            MGT_SET_CFA_3D
#define mgt_set_cfbx_3d           MGT_SET_CFBX_3D
#define mgt_set_cfby_3d           MGT_SET_CFBY_3D
#define mgt_set_cfbz_3d           MGT_SET_CFBZ_3D

#define mgt_set_cfs_2d            MGT_SET_CFS_2D
#define mgt_set_cfs_3d            MGT_SET_CFS_3D

#define mgt_get_vel_2d            MGT_GET_VEL_2D
#define mgt_set_vel_2d            MGT_SET_VEL_2D
#define mgt_get_vel_3d            MGT_GET_VEL_3D
#define mgt_set_vel_3d            MGT_SET_VEL_3D

#elif defined(BL_FORT_USE_UNDERSCORE)

#define mgt_init                  mgt_init_
#define mgt_alloc                 mgt_alloc_
#define mgt_nodal_alloc           mgt_nodal_alloc_
#define mgt_set_level             mgt_set_level_
#define mgt_set_nodal_level       mgt_set_nodal_level_
#define mgt_finalize              mgt_finalize_
#define mgt_nodal_finalize        mgt_nodal_finalize_
#define mgt_init_coeffs_lev       mgt_init_coeffs_lev_
#define mgt_init_nodal_coeffs_lev mgt_init_nodal_coeffs_lev_
#define mgt_finalize_stencil      mgt_finalize_stencil_
#define mgt_finalize_nodal_stencil mgt_finalize_nodal_stencil_
#define mgt_finalize_stencil_lev  mgt_finalize_stencil_lev_
#define mgt_finalize_nodal_stencil_lev  mgt_finalize_nodal_stencil_lev_
#define mgt_dealloc               mgt_dealloc_
#define mgt_nodal_dealloc         mgt_nodal_dealloc_
#define mgt_solve                 mgt_solve_
#define mgt_nodal_solve           mgt_nodal_solve_
#define mgt_divu                  mgt_divu_
#define mgt_newu                  mgt_newu_
#define mgt_set_defaults          mgt_set_defaults_
#define mgt_set_nodal_defaults    mgt_set_nodal_defaults_
#define mgt_get_defaults          mgt_get_defaults_
#define mgt_get_nodal_defaults    mgt_get_nodal_defaults_

#define mgt_set_rh_2d             mgt_set_rh_2d_
#define mgt_set_rh_3d             mgt_set_rh_3d_

#define mgt_set_uu_2d             mgt_set_uu_2d_
#define mgt_set_uu_3d             mgt_set_uu_3d_
#define mgt_get_uu_2d             mgt_get_uu_2d_
#define mgt_get_uu_3d             mgt_get_uu_3d_

#define mgt_set_pr_2d             mgt_set_pr_2d_
#define mgt_set_pr_3d             mgt_set_pr_3d_
#define mgt_get_pr_2d             mgt_get_pr_2d_
#define mgt_get_pr_3d             mgt_get_pr_3d_

#define mgt_set_cfa_2d            mgt_set_cfa_2d_
#define mgt_set_cfbx_2d           mgt_set_cfbx_2d_
#define mgt_set_cfby_2d           mgt_set_cfby_2d_
#define mgt_set_cfa_3d            mgt_set_cfa_3d_
#define mgt_set_cfbx_3d           mgt_set_cfbx_3d_
#define mgt_set_cfby_3d           mgt_set_cfby_3d_
#define mgt_set_cfbz_3d           mgt_set_cfbz_3d_

#define mgt_set_cfs_2d            mgt_set_cfs_2d_
#define mgt_set_cfs_3d            mgt_set_cfs_3d_

#define mgt_set_vel_2d            mgt_set_vel_2d_
#define mgt_get_vel_2d            mgt_get_vel_2d_
#define mgt_set_vel_3d            mgt_set_vel_3d_ 
#define mgt_get_vel_3d            mgt_get_vel_3d_

#elif defined(BL_FORT_USE_DBL_UNDERSCORE)
#define mgt_init                  mgt_init__
#define mgt_alloc                 mgt_alloc__
#define mgt_nodal_alloc           mgt_nodal_alloc__
#define mgt_set_level             mgt_set_level__
#define mgt_set_nodal_level       mgt_set_nodal_level__
#define mgt_finalize              mgt_finalize__
#define mgt_nodal_finalize        mgt_nodal_finalize__
#define mgt_init_coeffs_lev       mgt_init_coeffs_lev__
#define mgt_init_nodal_coeffs_lev mgt_init_nodal_coeffs_lev__
#define mgt_finalize_stencil      mgt_finalize_stencil__
#define mgt_finalize_nodal_stencil mgt_finalize_nodal_stencil__
#define mgt_finalize_stencil_lev  mgt_finalize_stencil_lev__
#define mgt_finalize_nodal_stencil_lev  mgt_finalize_nodal_stencil_lev__
#define mgt_dealloc               mgt_dealloc__
#define mgt_nodal_dealloc         mgt_nodal_dealloc__
#define mgt_solve                 mgt_solve__
#define mgt_divu                  mgt_divu__
#define mgt_newu                  mgt_newu__
#define mgt_nodal_solve           mgt_nodal_solve__
#define mgt_set_defaults          mgt_set_defaults__
#define mgt_set_nodal_defaults    mgt_set_nodal_defaults__
#define mgt_get_defaults          mgt_get_defaults__
#define mgt_get_nodal_defaults    mgt_get_nodal_defaults__

#define mgt_set_rh_2d             mgt_set_rh_2d__
#define mgt_set_rh_3d             mgt_set_rh_3d__

#define mgt_set_uu_2d             mgt_set_uu_2d__
#define mgt_get_uu_2d             mgt_get_uu_2d__
#define mgt_set_uu_3d             mgt_set_uu_3d__
#define mgt_get_uu_3d             mgt_get_uu_3d__

#define mgt_set_pr_2d             mgt_set_pr_2d__
#define mgt_get_pr_2d             mgt_get_pr_2d__
#define mgt_set_pr_3d             mgt_set_pr_3d__
#define mgt_get_pr_3d             mgt_get_pr_3d__

#define mgt_set_cfa_2d            mgt_set_cfa_2d__
#define mgt_set_cfbx_2d           mgt_set_cfbx_2d__
#define mgt_set_cfby_2d           mgt_set_cfby_2d__
#define mgt_set_cfa_3d            mgt_set_cfa_3d__
#define mgt_set_cfbx_3d           mgt_set_cfbx_3d__
#define mgt_set_cfby_3d           mgt_set_cfby_3d__
#define mgt_set_cfbz_3d           mgt_set_cfbz_3d__

#define mgt_set_cfs_2d            mgt_set_cfs_2d__
#define mgt_set_cfs_3d            mgt_set_cfs_3d__

#define mgt_set_vel_2d            mgt_set_vel_2d__
#define mgt_get_vel_2d            mgt_get_vel_2d__
#define mgt_set_vel_3d            mgt_set_vel_3d__
#define mgt_get_vel_3d            mgt_get_vel_3d__
#endif

#ifdef __cplusplus
extern "C" 
{
#endif

  /* The contants match the ones in 'bc.f90'; care is needed to ensure 
     that the continue to match. */

  const int MGT_BC_PER = -1;	/* Periodic  */
  const int MGT_BC_INT =  0;	/* Interior  */
  const int MGT_BC_DIR =  1;	/* Dirichlet */
  const int MGT_BC_NEU =  2;	/* Neumann   */

  void mgt_init();

  void mgt_alloc(const int* dm, const int* nlevel, const int* nodal);
  void mgt_nodal_alloc(const int* dm, const int* nlevel, const int* nodal, const int* stencil_type);

  void mgt_set_level(const int* lev, const int* nb, const int* dm, 
		     const int* lo, const int* hi, 
		     const int* pd_lo, const int* pd_hi, 
		     const int* bc, const int* pm,
		     const int* pmap);

  void mgt_set_nodal_level(const int* lev, const int* nb, const int* dm, 
		           const int* lo, const int* hi, 
      		           const int* pd_lo, const int* pd_hi, 
      		           const int* bc, const int* pm,
		           const int* pmap);

  void mgt_finalize(const Real* dx);
  void mgt_nodal_finalize(const Real* dx);

  void mgt_init_coeffs_lev(const int* lev);
  void mgt_init_nodal_coeffs_lev(const int* lev);

  void mgt_init_nodal_coeffs_lev(const int* lev);
  
  void mgt_finalize_stencil();
  void mgt_finalize_nodal_stencil();
  
  void mgt_finalize_stencil_lev(const int* lev,
			    const Real* xa, const Real* xb,
			    const Real* pxa, const Real* pxbb);

  void mgt_finalize_nodal_stencil_lev(const int* lev);
  
  void mgt_set_rh_2d(const int* lev, const int* n, const Real* rh, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_uu_2d(const int* lev, const int* n, Real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_set_uu_2d(const int* lev, const int* n, const Real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_pr_2d(const int* lev, const int* n, Real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_set_pr_2d(const int* lev, const int* n, const Real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_set_cfa_2d(const int* lev, const int* n, const Real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);
  
  void mgt_set_cfbx_2d(const int* lev, const int* n, const Real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  
  void mgt_set_cfby_2d(const int* lev, const int* n, const Real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  
  void mgt_set_cfs_2d(const int* lev, const int* n, const Real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);
  
  void mgt_get_vel_2d(const int* lev, const int* n, Real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);
  
  void mgt_set_vel_2d(const int* lev, const int* n, const Real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);



  
  void mgt_set_rh_3d(const int* lev, const int* n, const Real* rh, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_uu_3d(const int* lev, const int* n, Real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_set_uu_3d(const int* lev, const int* n, const Real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_pr_3d(const int* lev, const int* n, Real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_set_pr_3d(const int* lev, const int* n, const Real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_set_cfa_3d(const int* lev, const int* n, const Real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);
  
  void mgt_set_cfbx_3d(const int* lev, const int* n, const Real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  
  void mgt_set_cfby_3d(const int* lev, const int* n, const Real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  
  void mgt_set_cfbz_3d(const int* lev, const int* n, const Real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  
  void mgt_set_cfs_3d(const int* lev, const int* n, const Real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);
  
  void mgt_get_vel_3d(const int* lev, const int* n, Real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);
  
  void mgt_set_vel_3d(const int* lev, const int* n, const Real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);

  void mgt_dealloc();

  void mgt_nodal_dealloc();
  
  void mgt_solve(const Real& tol, const Real& abs_tol);
  
  void mgt_nodal_solve(const Real& tol, const Real& abs_tol);
  
  void mgt_divu();
  
  void mgt_newu();

  void mgt_set_defaults(const int* nu1, const int* nu2, const int* nub, const int* nuf,
                        const int* gamma, const Real* omega,
                        const int* max_iter, const int* bottom_max_iter,
                        const int* bottom_solver, const Real* bottom_solver_eps,
                        const int* verbose, const int* cg_verbose,
                        const int* max_nlevel, const int* min_width,
                        const int* cycle, const int* smoother, const int* stencil_type);

  void mgt_set_nodal_defaults(const int* nu1, const int* nu2, const int* nub, const int* nuf,
                              const int* gamma, const Real* omega,
                              const int* max_iter, const int* bottom_max_iter,
                              const int* bottom_solver, const Real* bottom_solver_eps,
                              const int* verbose, const int* cg_verbose,
                              const int* max_nlevel, const int* min_width,
                              const int* cycle, const int* smoother, const int* stencil_type);

  void mgt_get_defaults(const int* nu1, const int* nu2, const int* nub, const int* nuf,
                        const int* gamma, const Real* omega,
                        const int* max_iter, const int* bottom_max_iter,
                        const int* bottom_solver,
                        const int* verbose, const int* cg_verbose,
                        const int* max_nlevel, const int* min_width,
                        const int* cycle, const int* smoother);

  void mgt_get_nodal_defaults(const int* nu1, const int* nu2, const int* nub, const int* nuf,
                              const int* gamma, const Real* omega,
                              const int* max_iter, const int* bottom_max_iter,
                              const int* bottom_solver,
                              const int* verbose, const int* cg_verbose,
                              const int* max_nlevel, const int* min_width,
                              const int* cycle, const int* smoother);
  
#ifdef __cplusplus
}
#endif

#endif
