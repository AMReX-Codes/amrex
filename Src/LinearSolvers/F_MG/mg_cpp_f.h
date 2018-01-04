#ifndef _MG_CPP_F_H_
#define _MG_CPP_F_H_

#if defined(BL_FORT_USE_UPPERCASE)
#define mgt_init                  MGT_INIT
#define mgt_flush_copyassoc_cache MGT_FLUSH_COPYASSOC_CACHE
#define mgt_flush_output          MGT_FLUSH_OUTPUT
#define mgt_cc_alloc              MGT_CC_ALLOC
#define mgt_nodal_alloc           MGT_NODAL_ALLOC
#define mgt_set_level             MGT_SET_LEVEL
#define mgt_set_nodal_level       MGT_SET_NODAL_LEVEL
#define mgt_finalize              MGT_FINALIZE
#define mgt_finalize_n            MGT_FINALIZE_N
#define mgt_nodal_finalize        MGT_NODAL_FINALIZE
#define mgt_init_coeffs_lev       MGT_INIT_COEFFS_LEV
#define mgt_init_mc_coeffs_lev    MGT_INIT_MC_COEFFS_LEV
#define mgt_init_nodal_coeffs_lev MGT_INIT_NODAL_COEFFS_LEV
#define mgt_init_const_nodal_coeffs_lev MGT_INIT_CONST_NODAL_COEFFS_LEV
#define mgt_finalize_stencil      MGT_FINALIZE_STENCIL
#define mgt_finalize_nodal_stencil MGT_FINALIZE_NODAL_STENCIL
#define mgt_finalize_stencil_lev  MGT_FINALIZE_STENCIL_LEV
#define mgt_finalize_const_stencil_lev  MGT_FINALIZE_CONST_STENCIL_LEV
#define mgt_mc_finalize_stencil_lev  MGT_MC_FINALIZE_STENCIL_LEV
#define mgt_finalize_nodal_stencil_lev  MGT_FINALIZE_NODAL_STENCIL_LEV
#define mgt_dealloc               MGT_DEALLOC
#define mgt_nodal_dealloc         MGT_NODAL_DEALLOC
#define mgt_applyop               MGT_APPLYOP
#define mgt_solve                 MGT_SOLVE
#define mgt_compute_flux          MGT_COMPUTE_FLUX
#define mgt_delete_flux           MGT_DELETE_FLUX
#define mgt_compute_residual      MGT_COMPUTE_RESIDUAL
#define mgt_nodal_solve           MGT_NODAL_SOLVE
#define mgt_divu                  MGT_DIVU
#define mgt_newu                  MGT_NEWU
#define mgt_set_defaults          MGT_SET_DEFAULTS
#define mgt_set_nodal_defaults    MGT_SET_NODAL_DEFAULTS
#define mgt_get_defaults          MGT_GET_DEFAULTS
#define mgt_get_nodal_defaults    MGT_GET_NODAL_DEFAULTS
#define mgt_set_maxorder          MGT_SET_MAXORDER

#define mgt_set_rh_1d             MGT_SET_RH_1D
#define mgt_set_rh_2d             MGT_SET_RH_2D
#define mgt_set_rh_3d             MGT_SET_RH_3D

#define mgt_add_rh_nodal_1d       MGT_ADD_RH_NODAL_1D
#define mgt_add_rh_nodal_2d       MGT_ADD_RH_NODAL_2D
#define mgt_add_rh_nodal_3d       MGT_ADD_RH_NODAL_3D

#define mgt_set_rh_nodal_1d       MGT_SET_RH_NODAL_1D
#define mgt_set_rh_nodal_2d       MGT_SET_RH_NODAL_2D
#define mgt_set_rh_nodal_3d       MGT_SET_RH_NODAL_3D

#define mgt_set_uu_1d             MGT_SET_UU_1D
#define mgt_get_uu_1d             MGT_GET_UU_1D
#define mgt_set_uu_2d             MGT_SET_UU_2D
#define mgt_get_uu_2d             MGT_GET_UU_2D
#define mgt_set_uu_3d             MGT_SET_UU_3D
#define mgt_get_uu_3d             MGT_GET_UU_3D

#define mgt_get_res_1d            MGT_GET_RES_1D
#define mgt_get_res_2d            MGT_GET_RES_2D
#define mgt_get_res_3d            MGT_GET_RES_3D

#define mgt_set_pr_1d             MGT_SET_PR_1D
#define mgt_get_pr_1d             MGT_GET_PR_1D
#define mgt_set_pr_2d             MGT_SET_PR_2D
#define mgt_get_pr_2d             MGT_GET_PR_2D
#define mgt_set_pr_3d             MGT_SET_PR_3D
#define mgt_get_pr_3d             MGT_GET_PR_3D

#define mgt_get_gp_1d             MGT_SET_GP_1D
#define mgt_get_gp_2d             MGT_SET_GP_2D
#define mgt_get_gp_3d             MGT_SET_GP_3D

#define mgt_set_cfa_1d            MGT_SET_CFA_1D
#define mgt_set_cfaa_1d           MGT_SET_CFAA_1D
#define mgt_set_cfa2_1d           MGT_SET_CFA2_1D
#define mgt_set_cfbx_1d           MGT_SET_CFBX_1D
#define mgt_set_cfbnx_1d          MGT_SET_CFBNX_1D
#define mgt_set_cfa_2d            MGT_SET_CFA_2D
#define mgt_set_cfaa_2d           MGT_SET_CFAA_2D
#define mgt_set_cfa2_2d           MGT_SET_CFA2_2D
#define mgt_set_cfbx_2d           MGT_SET_CFBX_2D
#define mgt_set_cfby_2d           MGT_SET_CFBY_2D
#define mgt_set_cfbnx_2d          MGT_SET_CFBNX_2D
#define mgt_set_cfbny_2d          MGT_SET_CFBNY_2D
#define mgt_set_cfa2_3d           MGT_SET_CFA2_3D
#define mgt_set_cfbx_3d           MGT_SET_CFBX_3D
#define mgt_set_cfby_3d           MGT_SET_CFBY_3D
#define mgt_set_cfbz_3d           MGT_SET_CFBZ_3D
#define mgt_set_cfbnx_3d          MGT_SET_CFBNX_3D
#define mgt_set_cfbny_3d          MGT_SET_CFBNY_3D
#define mgt_set_cfbnz_3d          MGT_SET_CFBNZ_3D

#define mgt_set_cfa_1d_const      MGT_SET_CFA_1D_CONST
#define mgt_set_cfbx_1d_const     MGT_SET_CFBX_1D_CONST
#define mgt_set_cfa_2d_const      MGT_SET_CFA_2D_CONST
#define mgt_set_cfbx_2d_const     MGT_SET_CFBX_2D_CONST
#define mgt_set_cfby_2d_const     MGT_SET_CFBY_2D_CONST
#define mgt_set_cfa_3d_const      MGT_SET_CFA_3D_CONST
#define mgt_set_cfbx_3d_const     MGT_SET_CFBX_3D_CONST
#define mgt_set_cfby_3d_const     MGT_SET_CFBY_3D_CONST
#define mgt_set_cfbz_3d_const     MGT_SET_CFBZ_3D_CONST

#define mgt_set_cfs_1d            MGT_SET_CFS_1D
#define mgt_set_cfs_2d            MGT_SET_CFS_2D
#define mgt_set_cfs_3d            MGT_SET_CFS_3D

#define mgt_set_cfs_1d_const      MGT_SET_CFS_1D_CONST
#define mgt_set_cfs_2d_const      MGT_SET_CFS_2D_CONST
#define mgt_set_cfs_3d_const      MGT_SET_CFS_3D_CONST

#define mgt_get_vel_1d            MGT_GET_VEL_1D
#define mgt_set_vel_1d            MGT_SET_VEL_1D
#define mgt_get_vel_2d            MGT_GET_VEL_2D
#define mgt_set_vel_2d            MGT_SET_VEL_2D
#define mgt_get_vel_3d            MGT_GET_VEL_3D
#define mgt_set_vel_3d            MGT_SET_VEL_3D

#define mgt_alloc_nodal_sync      MGT_ALLOC_NODAL_SYNC
#define mgt_dealloc_nodal_sync    MGT_DEALLOC_NODAL_SYNC
#define mgt_compute_sync_resid_crse  MGT_COMPUTE_SYNC_RESID_CRSE
#define mgt_compute_sync_resid_fine  MGT_COMPUTE_SYNC_RESID_FINE

#define mgt_set_sync_msk_1d       MGT_SET_SYNC_MSK_1D 
#define mgt_set_sync_msk_2d       MGT_SET_SYNC_MSK_2D 
#define mgt_set_sync_msk_3d       MGT_SET_SYNC_MSK_3D 

#define mgt_set_vold_1d           MGT_SET_VOLD_1D 
#define mgt_set_vold_2d           MGT_SET_VOLD_2D 
#define mgt_set_vold_3d           MGT_SET_VOLD_3D 

#define mgt_get_sync_res_1d       MGT_GET_SYNC_RES_1D 
#define mgt_get_sync_res_2d       MGT_GET_SYNC_RES_2D 
#define mgt_get_sync_res_3d       MGT_GET_SYNC_RES_3D 

#define mgt_alloc_rhcc_nodal      MGT_ALLOC_RHCC_NODAL
#define mgt_dealloc_rhcc_nodal    MGT_DEALLOC_RHCC_NODAL

#define mgt_set_rhcc_nodal_1d     MGT_SET_RHCC_NODAL_1D
#define mgt_set_rhcc_nodal_2d     MGT_SET_RHCC_NODAL_2D
#define mgt_set_rhcc_nodal_3d     MGT_SET_RHCC_NODAL_3D

#define mgt_add_divucc            MGT_ADD_DIVUCC

#elif defined(BL_FORT_USE_UNDERSCORE)

#define mgt_init                  mgt_init_
#define mgt_flush_copyassoc_cache mgt_flush_copyassoc_cache_
#define mgt_flush_output          mgt_flush_output_
#define mgt_alloc                 mgt_alloc_
#define mgt_cc_alloc              mgt_cc_alloc_
#define mgt_nodal_alloc           mgt_nodal_alloc_
#define mgt_set_level             mgt_set_level_
#define mgt_set_nodal_level       mgt_set_nodal_level_
#define mgt_finalize              mgt_finalize_
#define mgt_finalize_n            mgt_finalize_n_
#define mgt_nodal_finalize        mgt_nodal_finalize_
#define mgt_init_coeffs_lev       mgt_init_coeffs_lev_
#define mgt_init_mc_coeffs_lev    mgt_init_mc_coeffs_lev_
#define mgt_init_nodal_coeffs_lev mgt_init_nodal_coeffs_lev_
#define mgt_init_const_nodal_coeffs_lev mgt_init_const_nodal_coeffs_lev_
#define mgt_finalize_stencil      mgt_finalize_stencil_
#define mgt_finalize_nodal_stencil mgt_finalize_nodal_stencil_
#define mgt_finalize_stencil_lev  mgt_finalize_stencil_lev_
#define mgt_finalize_const_stencil_lev  mgt_finalize_const_stencil_lev_
#define mgt_mc_finalize_stencil_lev   mgt_mc_finalize_stencil_lev_
#define mgt_finalize_nodal_stencil_lev  mgt_finalize_nodal_stencil_lev_
#define mgt_dealloc               mgt_dealloc_
#define mgt_nodal_dealloc         mgt_nodal_dealloc_
#define mgt_solve                 mgt_solve_
#define mgt_applyop               mgt_applyop_
#define mgt_compute_flux          mgt_compute_flux_
#define mgt_delete_flux           mgt_delete_flux_
#define mgt_compute_residual      mgt_compute_residual_
#define mgt_nodal_solve           mgt_nodal_solve_
#define mgt_divu                  mgt_divu_
#define mgt_newu                  mgt_newu_
#define mgt_set_defaults          mgt_set_defaults_
#define mgt_set_nodal_defaults    mgt_set_nodal_defaults_
#define mgt_get_defaults          mgt_get_defaults_
#define mgt_get_nodal_defaults    mgt_get_nodal_defaults_
#define mgt_set_maxorder          mgt_set_maxorder_

#define mgt_set_rh_1d             mgt_set_rh_1d_
#define mgt_set_rh_2d             mgt_set_rh_2d_
#define mgt_set_rh_3d             mgt_set_rh_3d_

#define mgt_add_rh_nodal_1d       mgt_add_rh_nodal_1d_
#define mgt_add_rh_nodal_2d       mgt_add_rh_nodal_2d_
#define mgt_add_rh_nodal_3d       mgt_add_rh_nodal_3d_

#define mgt_set_rh_nodal_1d       mgt_set_rh_nodal_1d_
#define mgt_set_rh_nodal_2d       mgt_set_rh_nodal_2d_
#define mgt_set_rh_nodal_3d       mgt_set_rh_nodal_3d_

#define mgt_set_uu_1d             mgt_set_uu_1d_
#define mgt_get_uu_1d             mgt_get_uu_1d_
#define mgt_set_uu_2d             mgt_set_uu_2d_
#define mgt_get_uu_2d             mgt_get_uu_2d_
#define mgt_set_uu_3d             mgt_set_uu_3d_
#define mgt_get_uu_3d             mgt_get_uu_3d_

#define mgt_get_res_1d            mgt_get_res_1d_
#define mgt_get_res_2d            mgt_get_res_2d_
#define mgt_get_res_3d            mgt_get_res_3d_

#define mgt_set_pr_1d             mgt_set_pr_1d_
#define mgt_set_pr_2d             mgt_set_pr_2d_
#define mgt_set_pr_3d             mgt_set_pr_3d_
#define mgt_get_pr_1d             mgt_get_pr_1d_
#define mgt_get_pr_2d             mgt_get_pr_2d_
#define mgt_get_pr_3d             mgt_get_pr_3d_

#define mgt_get_gp_1d             mgt_get_gp_1d_
#define mgt_get_gp_2d             mgt_get_gp_2d_
#define mgt_get_gp_3d             mgt_get_gp_3d_

#define mgt_set_cfa_1d            mgt_set_cfa_1d_
#define mgt_set_cfaa_1d           mgt_set_cfaa_1d_
#define mgt_set_cfa2_1d           mgt_set_cfa2_1d_
#define mgt_set_cfbx_1d           mgt_set_cfbx_1d_
#define mgt_set_cfbnx_1d          mgt_set_cfbnx_1d_
#define mgt_set_cfnbx_1d          mgt_set_cfbnx_1d_
#define mgt_set_cfa_2d            mgt_set_cfa_2d_
#define mgt_set_cfaa_2d           mgt_set_cfaa_2d_
#define mgt_set_cfa2_2d           mgt_set_cfa2_2d_
#define mgt_set_cfbx_2d           mgt_set_cfbx_2d_
#define mgt_set_cfby_2d           mgt_set_cfby_2d_
#define mgt_set_cfbnx_2d          mgt_set_cfbnx_2d_
#define mgt_set_cfbny_2d          mgt_set_cfbny_2d_
#define mgt_set_cfa_3d            mgt_set_cfa_3d_
#define mgt_set_cfaa_3d           mgt_set_cfaa_3d_
#define mgt_set_cfa2_3d           mgt_set_cfa2_3d_
#define mgt_set_cfbx_3d           mgt_set_cfbx_3d_
#define mgt_set_cfby_3d           mgt_set_cfby_3d_
#define mgt_set_cfbz_3d           mgt_set_cfbz_3d_
#define mgt_set_cfbnx_3d          mgt_set_cfbnx_3d_
#define mgt_set_cfbny_3d          mgt_set_cfbny_3d_
#define mgt_set_cfbnz_3d          mgt_set_cfbnz_3d_

#define mgt_set_cfa_1d_const      mgt_set_cfa_1d_const_
#define mgt_set_cfbx_1d_const     mgt_set_cfbx_1d_const_
#define mgt_set_cfa_2d_const      mgt_set_cfa_2d_const_
#define mgt_set_cfbx_2d_const     mgt_set_cfbx_2d_const_
#define mgt_set_cfby_2d_const     mgt_set_cfby_2d_const_
#define mgt_set_cfa_3d_const      mgt_set_cfa_3d_const_
#define mgt_set_cfbx_3d_const     mgt_set_cfbx_3d_const_
#define mgt_set_cfby_3d_const     mgt_set_cfby_3d_const_
#define mgt_set_cfbz_3d_const     mgt_set_cfbz_3d_const_

#define mgt_set_cfs_1d            mgt_set_cfs_1d_
#define mgt_set_cfs_2d            mgt_set_cfs_2d_
#define mgt_set_cfs_3d            mgt_set_cfs_3d_

#define mgt_set_vel_1d            mgt_set_vel_1d_
#define mgt_get_vel_1d            mgt_get_vel_1d_
#define mgt_set_vel_2d            mgt_set_vel_2d_
#define mgt_get_vel_2d            mgt_get_vel_2d_
#define mgt_set_vel_3d            mgt_set_vel_3d_ 
#define mgt_get_vel_3d            mgt_get_vel_3d_

#define mgt_alloc_nodal_sync      mgt_alloc_nodal_sync_
#define mgt_dealloc_nodal_sync    mgt_dealloc_nodal_sync_
#define mgt_compute_sync_resid_crse  mgt_compute_sync_resid_crse_
#define mgt_compute_sync_resid_fine  mgt_compute_sync_resid_fine_

#define mgt_set_sync_msk_1d       mgt_set_sync_msk_1d_ 
#define mgt_set_sync_msk_2d       mgt_set_sync_msk_2d_ 
#define mgt_set_sync_msk_3d       mgt_set_sync_msk_3d_ 

#define mgt_set_vold_1d           mgt_set_vold_1d_
#define mgt_set_vold_2d           mgt_set_vold_2d_
#define mgt_set_vold_3d           mgt_set_vold_3d_

#define mgt_get_sync_res_1d       mgt_get_sync_res_1d_ 
#define mgt_get_sync_res_2d       mgt_get_sync_res_2d_ 
#define mgt_get_sync_res_3d       mgt_get_sync_res_3d_ 

#define mgt_alloc_rhcc_nodal      mgt_alloc_rhcc_nodal_
#define mgt_dealloc_rhcc_nodal    mgt_dealloc_rhcc_nodal_

#define mgt_set_rhcc_nodal_1d     mgt_set_rhcc_nodal_1d_
#define mgt_set_rhcc_nodal_2d     mgt_set_rhcc_nodal_2d_
#define mgt_set_rhcc_nodal_3d     mgt_set_rhcc_nodal_3d_

#define mgt_add_divucc            mgt_add_divucc_

#elif defined(BL_FORT_USE_DBL_UNDERSCORE)

#define mgt_init                  mgt_init__
#define mgt_flush_copyassoc_cache mgt_flush_copyassoc_cache__
#define mgt_flush_output          mgt_flush_output__
#define mgt_cc_alloc              mgt_cc_alloc_
#define mgt_nodal_alloc           mgt_nodal_alloc__
#define mgt_set_level             mgt_set_level__
#define mgt_set_nodal_level       mgt_set_nodal_level__
#define mgt_finalize              mgt_finalize__
#define mgt_finalize_n            mgt_finalize_n__
#define mgt_nodal_finalize        mgt_nodal_finalize__
#define mgt_init_coeffs_lev       mgt_init_coeffs_lev__
#define mgt_init_mc_coeffs_lev    mgt_init_mc_coeffs_lev__
#define mgt_init_nodal_coeffs_lev mgt_init_nodal_coeffs_lev__
#define mgt_init_const_nodal_coeffs_lev mgt_init_const_nodal_coeffs_lev__
#define mgt_finalize_stencil      mgt_finalize_stencil__
#define mgt_finalize_nodal_stencil mgt_finalize_nodal_stencil__
#define mgt_finalize_stencil_lev  mgt_finalize_stencil_lev__
#define mgt_finalize_const_stencil_lev  mgt_finalize_const_stencil_lev__
#define mgt_mc_finalize_stencil_lev   mgt_mc_finalize_stencil_lev__
#define mgt_finalize_nodal_stencil_lev  mgt_finalize_nodal_stencil_lev__
#define mgt_dealloc               mgt_dealloc__
#define mgt_nodal_dealloc         mgt_nodal_dealloc__
#define mgt_solve                 mgt_solve__
#define mgt_applyop               mgt_applyop__
#define mgt_compute_flux          mgt_compute_flux__
#define mgt_delete_flux           mgt_delete_flux_
#define mgt_compute_residual      mgt_compute_residual__
#define mgt_divu                  mgt_divu__
#define mgt_newu                  mgt_newu__
#define mgt_nodal_solve           mgt_nodal_solve__
#define mgt_set_defaults          mgt_set_defaults__
#define mgt_set_nodal_defaults    mgt_set_nodal_defaults__
#define mgt_get_defaults          mgt_get_defaults__
#define mgt_get_nodal_defaults    mgt_get_nodal_defaults__
#define mgt_set_maxorder          mgt_set_maxorder_

#define mgt_set_rh_1d             mgt_set_rh_1d__
#define mgt_set_rh_2d             mgt_set_rh_2d__
#define mgt_set_rh_3d             mgt_set_rh_3d__

#define mgt_add_rh_nodal_1d       mgt_add_rh_nodal_1d__
#define mgt_add_rh_nodal_2d       mgt_add_rh_nodal_2d__
#define mgt_add_rh_nodal_3d       mgt_add_rh_nodal_3d__

#define mgt_set_rh_nodal_1d       mgt_set_rh_nodal_1d__
#define mgt_set_rh_nodal_2d       mgt_set_rh_nodal_2d__
#define mgt_set_rh_nodal_3d       mgt_set_rh_nodal_3d__

#define mgt_set_uu_1d             mgt_set_uu_1d__
#define mgt_get_uu_1d             mgt_get_uu_1d__
#define mgt_set_uu_2d             mgt_set_uu_2d__
#define mgt_get_uu_2d             mgt_get_uu_2d__
#define mgt_set_uu_3d             mgt_set_uu_3d__
#define mgt_get_uu_3d             mgt_get_uu_3d__

#define mgt_get_res_1d            mgt_get_res_1d__
#define mgt_get_res_2d            mgt_get_res_2d__
#define mgt_get_res_3d            mgt_get_res_3d__

#define mgt_set_pr_1d             mgt_set_pr_1d__
#define mgt_get_pr_1d             mgt_get_pr_1d__
#define mgt_set_pr_2d             mgt_set_pr_2d__
#define mgt_get_pr_2d             mgt_get_pr_2d__
#define mgt_set_pr_3d             mgt_set_pr_3d__
#define mgt_get_pr_3d             mgt_get_pr_3d__

#define mgt_get_gp_1d             mgt_get_gp_1d__
#define mgt_get_gp_2d             mgt_get_gp_2d__
#define mgt_set_gp_3d             mgt_get_gp_3d__

#define mgt_set_cfa_1d            mgt_set_cfa_1d__
#define mgt_set_cfaa_1d           mgt_set_cfaa_1d__
#define mgt_set_cfa2_1d           mgt_set_cfa2_1d__
#define mgt_set_cfbx_1d           mgt_set_cfbx_1d__
#define mgt_set_cfa_2d            mgt_set_cfa_2d__
#define mgt_set_cfaa_2d           mgt_set_cfaa_2d__
#define mgt_set_cfa2_2d           mgt_set_cfa2_2d__
#define mgt_set_cfbx_2d           mgt_set_cfbx_2d__
#define mgt_set_cfby_2d           mgt_set_cfby_2d__
#define mgt_set_cfa_3d            mgt_set_cfa_3d__
#define mgt_set_cfaa_3d           mgt_set_cfaa_3d__
#define mgt_set_cfa2_3d           mgt_set_cfa2_3d__
#define mgt_set_cfbx_3d           mgt_set_cfbx_3d__
#define mgt_set_cfby_3d           mgt_set_cfby_3d__
#define mgt_set_cfbz_3d           mgt_set_cfbz_3d__

#define mgt_set_cfa_1d_const      mgt_set_cfa_1d_const__
#define mgt_set_cfbx_1d_const     mgt_set_cfbx_1d_const__
#define mgt_set_cfa_2d_const      mgt_set_cfa_2d_const__
#define mgt_set_cfbx_2d_const     mgt_set_cfbx_2d_const__
#define mgt_set_cfby_2d_const     mgt_set_cfby_2d_const__
#define mgt_set_cfa_3d_const      mgt_set_cfa_3d_const__
#define mgt_set_cfbx_3d_const     mgt_set_cfbx_3d_const__
#define mgt_set_cfby_3d_const     mgt_set_cfby_3d_const__
#define mgt_set_cfbz_3d_const     mgt_set_cfbz_3d_const__

#define mgt_set_cfs_1d            mgt_set_cfs_1d__
#define mgt_set_cfs_2d            mgt_set_cfs_2d__
#define mgt_set_cfs_3d            mgt_set_cfs_3d__

#define mgt_set_vel_1d            mgt_set_vel_2d__
#define mgt_get_vel_2d            mgt_get_vel_1d__
#define mgt_set_vel_2d            mgt_set_vel_2d__
#define mgt_get_vel_2d            mgt_get_vel_2d__
#define mgt_set_vel_3d            mgt_set_vel_3d__
#define mgt_get_vel_3d            mgt_get_vel_3d__

#define mgt_alloc_nodal_sync      mgt_alloc_nodal_sync__
#define mgt_dealloc_nodal_sync    mgt_dealloc_nodal_sync__
#define mgt_compute_sync_resid_crse  mgt_compute_sync_resid_crse__
#define mgt_compute_sync_resid_fine  mgt_compute_sync_resid_fine__

#define mgt_set_sync_msk_1d       mgt_set_sync_msk_1d__ 
#define mgt_set_sync_msk_2d       mgt_set_sync_msk_2d__ 
#define mgt_set_sync_msk_3d       mgt_set_sync_msk_3d__ 

#define mgt_set_vold_1d           mgt_set_vold_1d__
#define mgt_set_vold_2d           mgt_set_vold_2d__
#define mgt_set_vold_3d           mgt_set_vold_3d__

#define mgt_get_sync_res_1d       mgt_get_sync_res_1d__ 
#define mgt_get_sync_res_2d       mgt_get_sync_res_2d__ 
#define mgt_get_sync_res_3d       mgt_get_sync_res_3d__ 

#define mgt_alloc_rhcc_nodal      mgt_alloc_rhcc_nodal__
#define mgt_dealloc_rhcc_nodal    mgt_dealloc_rhcc_nodal__

#define mgt_set_rhcc_nodal_1d     mgt_set_rhcc_nodal_1d_
#define mgt_set_rhcc_nodal_2d     mgt_set_rhcc_nodal_2d_
#define mgt_set_rhcc_nodal_3d     mgt_set_rhcc_nodal_3d_

#define mgt_add_divucc            mgt_add_divucc_

#endif

#ifdef __cplusplus
extern "C" 
{
#endif

  /* The contants match the ones in 'bc.f90'; care is needed to ensure 
     that they continue to match. */

  const int MGT_BC_PER = -1;	/* Periodic  */
  const int MGT_BC_INT =  0;	/* Interior  */
  const int MGT_BC_DIR =  1;	/* Dirichlet */
  const int MGT_BC_NEU =  2;	/* Neumann   */

  void mgt_init();

  void mgt_flush_copyassoc_cache();

  void mgt_flush_output();

  void mgt_cc_alloc   (const int* dm, const int* nlevel, const int* stencil_type);
  void mgt_nodal_alloc(const int* dm, const int* nlevel, const int* stencil_type);

  void mgt_set_level(const int* lev, const int* nb, const int* dm, 
		     const int* lo, const int* hi, 
		     const int* pd_lo, const int* pd_hi, 
		     const int* pm,
		     const int* pmap);

  void mgt_set_nodal_level(const int* lev, const int* nb, const int* dm, 
		           const int* lo, const int* hi, 
      		           const int* pd_lo, const int* pd_hi, 
      		           const int* pm,
		           const int* pmap);

  void mgt_finalize(const amrex_real* dx, const int* bc);
  void mgt_finalize_n(const amrex_real* dx, const int* bc, const int* nc_in);
  void mgt_nodal_finalize(const amrex_real* dx, const int* bc);

  void mgt_init_coeffs_lev(const int* lev);
  void mgt_init_mc_coeffs_lev(const int* lev, const int* nc, int* nc_opt);
  void mgt_init_nodal_coeffs_lev(const int* lev);
  void mgt_init_const_nodal_coeffs_lev(const int* lev, const amrex_real* val);
  
  void mgt_finalize_stencil();
  void mgt_finalize_nodal_stencil();
  
  void mgt_finalize_stencil_lev(const int* lev,
			    const amrex_real* xa, const amrex_real* xb,
			    const amrex_real* pxa, const amrex_real* pxb,
			    const int* dm);
  
  void mgt_finalize_const_stencil_lev(const int* lev,
			         const amrex_real* alpha_const, const amrex_real* beta_const,
			         const amrex_real* xa, const amrex_real* xb,
     			         const amrex_real* pxa, const amrex_real* pxb,
         			 const int* dm);

  void mgt_mc_finalize_stencil_lev(const int* lev,
				   const amrex_real* xa, const amrex_real* xb,
				   const amrex_real* pxa, const amrex_real* pxb,
				   const int* dm,
			    	   const int* nc_opt);

  void mgt_finalize_nodal_stencil_lev(const int* lev);
    
  void mgt_set_rh_1d(const int* lev, const int* n, const amrex_real* rh, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_add_rh_nodal_1d(const int* lev, const int* n, const amrex_real* rh, 
                           const int* plo, const int* phi, 
                           const int* lo, const int* hi,
                           const amrex_real* rhmax);

  void mgt_set_rh_nodal_1d(const int* lev, const int* n, const amrex_real* rh, 
                           const int* plo, const int* phi, 
                           const int* lo, const int* hi);
  
  void mgt_get_uu_1d(const int* lev, const int* n, amrex_real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_res_1d(const int* lev, const int* n, amrex_real* uu, 
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);
  
  void mgt_set_uu_1d(const int* lev, const int* n, const amrex_real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_gp_1d(const int* lev, const int* dir, const int* n,
		     amrex_real* gp,
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_pr_1d(const int* lev, const int* n, amrex_real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi,
		     const int& np, const int& ip);
  
  void mgt_set_pr_1d(const int* lev, const int* n, const amrex_real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi,
		     const int& np, const int& ip);
  
  void mgt_set_cfa_1d(const int* lev, const int* n, const amrex_real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);

  void mgt_set_cfaa_1d(const int* lev, const int* n, const amrex_real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi, const amrex_real* a);
  
  void mgt_set_cfa2_1d(const int* lev, const int* n, const amrex_real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi, const int& ncomps);
  
  void mgt_set_cfbx_1d(const int* lev, const int* n, const amrex_real* cf,
		       const amrex_real* b,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  
  void mgt_set_cfbnx_1d(const int* lev, const int* n, const amrex_real* cf,
		        const amrex_real* b,
		        const int* plo, const int* phi, 
		        const int* lo, const int* hi, const int& ncomps);
  
  void mgt_set_cfa_1d_const(const int* lev, const int* n, 
		            const int* lo, const int* hi, const amrex_real* value);
  
  void mgt_set_cfbx_1d_const(const int* lev, const int* n, 
		             const int* lo, const int* hi, const amrex_real* value);
  
  void mgt_set_cfs_1d(const int* lev, const int* n, const amrex_real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);
  
  void mgt_get_vel_1d(const int* lev, const int* n, amrex_real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi,
		      const int& nv, const int& iv);
  
  void mgt_set_vel_1d(const int* lev, const int* n, const amrex_real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi,
		      const int& nv, const int& iv);

  void mgt_set_rh_2d(const int* lev, const int* n, const amrex_real* rh, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);

  void mgt_add_rh_nodal_2d(const int* lev, const int* n, const amrex_real* rh, 
                           const int* plo, const int* phi, 
                           const int* lo, const int* hi,
                           const amrex_real* rhmax);
  
  void mgt_set_rh_nodal_2d(const int* lev, const int* n, const amrex_real* rh, 
                           const int* plo, const int* phi, 
                           const int* lo, const int* hi);

  void mgt_get_uu_2d(const int* lev, const int* n, amrex_real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_res_2d(const int* lev, const int* n, amrex_real* uu, 
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);
  
  void mgt_set_uu_2d(const int* lev, const int* n, const amrex_real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_gp_2d(const int* lev, const int* dir, const int* n,
		     amrex_real* gp,
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_pr_2d(const int* lev, const int* n, amrex_real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi,
		     const int& np, const int& ip);
  
  void mgt_set_pr_2d(const int* lev, const int* n, const amrex_real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi,
		     const int& np, const int& ip);
  
  void mgt_set_cfa_2d(const int* lev, const int* n, const amrex_real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);

  void mgt_set_cfaa_2d(const int* lev, const int* n, const amrex_real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi, const amrex_real* a);
  
  void mgt_set_cfa2_2d(const int* lev, const int* n, const amrex_real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi, const int& ncomps);
  
  void mgt_set_cfbx_2d(const int* lev, const int* n, const amrex_real* cf,
		       const amrex_real* b,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  
  void mgt_set_cfby_2d(const int* lev, const int* n, const amrex_real* cf,
		       const amrex_real* b,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  
  void mgt_set_cfbnx_2d(const int* lev, const int* n, const amrex_real* cf,
		        const amrex_real* b,
		        const int* plo, const int* phi, 
		        const int* lo, const int* hi, const int& ncomps);
  
  void mgt_set_cfbny_2d(const int* lev, const int* n, const amrex_real* cf,
		        const amrex_real* b,
		        const int* plo, const int* phi, 
		        const int* lo, const int* hi, const int& ncomps);
  
  void mgt_set_cfa_2d_const(const int* lev, const int* n, 
		            const int* lo, const int* hi, const amrex_real* value);
  
  void mgt_set_cfbx_2d_const(const int* lev, const int* n, 
		             const int* lo, const int* hi, const amrex_real* value);
  
  void mgt_set_cfby_2d_const(const int* lev, const int* n,
		             const int* lo, const int* hi, const amrex_real* value);
  
  void mgt_set_cfs_2d(const int* lev, const int* n, const amrex_real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);
  
  void mgt_get_vel_2d(const int* lev, const int* n, amrex_real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi,
		      const int& nv, const int& iv);
  
  void mgt_set_vel_2d(const int* lev, const int* n, const amrex_real* v,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi,
		      const int& nv, const int& iv);

  void mgt_set_rh_3d(const int* lev, const int* n, const amrex_real* rh, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);

  void mgt_add_rh_nodal_3d(const int* lev, const int* n, const amrex_real* rh, 
                           const int* plo, const int* phi, 
                           const int* lo, const int* hi,
                           const amrex_real* rhmax);
  
  void mgt_set_rh_nodal_3d(const int* lev, const int* n, const amrex_real* rh, 
                           const int* plo, const int* phi, 
                           const int* lo, const int* hi);

  void mgt_get_uu_3d(const int* lev, const int* n, amrex_real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_set_uu_3d(const int* lev, const int* n, const amrex_real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_res_3d(const int* lev, const int* n, amrex_real* uu, 
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);
  
  void mgt_get_gp_3d(const int* lev, const int* dir, const int* n, 
		     amrex_real* gp,
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_pr_3d(const int* lev, const int* n, amrex_real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi,
		     const int& np, const int& ip);
  
  void mgt_set_pr_3d(const int* lev, const int* n, const amrex_real* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi,
		     const int& np, const int& ip);
  
  void mgt_set_cfa_3d(const int* lev, const int* n, const amrex_real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);

  void mgt_set_cfaa_3d(const int* lev, const int* n, const amrex_real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi, const amrex_real* a);
  
  void mgt_set_cfa2_3d(const int* lev, const int* n, const amrex_real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi, const int& ncomps);
  
  void mgt_set_cfbx_3d(const int* lev, const int* n, const amrex_real* cf,
		       const amrex_real* b,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  
  void mgt_set_cfby_3d(const int* lev, const int* n, const amrex_real* cf,
		       const amrex_real* b,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  
  void mgt_set_cfbz_3d(const int* lev, const int* n, const amrex_real* cf,
		       const amrex_real* b,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  
  void mgt_set_cfbnx_3d(const int* lev, const int* n, const amrex_real* cf,
		        const amrex_real* b,
		        const int* plo, const int* phi, 
		        const int* lo, const int* hi, const int& ncomps);
  
  void mgt_set_cfbny_3d(const int* lev, const int* n, const amrex_real* cf,
		        const amrex_real* b,
		        const int* plo, const int* phi, 
		        const int* lo, const int* hi, const int& ncomps);
  
  void mgt_set_cfbnz_3d(const int* lev, const int* n, const amrex_real* cf,
		        const amrex_real* b,
		        const int* plo, const int* phi, 
		        const int* lo, const int* hi, const int& ncomps);
  
  void mgt_set_cfa_3d_const(const int* lev, const int* n, 
		            const int* lo, const int* hi, const amrex_real* value);
  
  void mgt_set_cfbx_3d_const(const int* lev, const int* n, 
		             const int* lo, const int* hi, const amrex_real* value);
  
  void mgt_set_cfby_3d_const(const int* lev, const int* n, 
		             const int* lo, const int* hi, const amrex_real* value);
  
  void mgt_set_cfbz_3d_const(const int* lev, const int* n, 
		             const int* lo, const int* hi, const amrex_real* value);
  
  void mgt_set_cfs_3d(const int* lev, const int* n, const amrex_real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi);
  
  void mgt_get_vel_3d(const int* lev, const int* n, amrex_real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi,
		      const int& nv, const int& iv);
  
  void mgt_set_vel_3d(const int* lev, const int* n, const amrex_real* cf,
		      const int* plo, const int* phi, 
		      const int* lo, const int* hi,
		      const int& nv, const int& iv);

  void mgt_set_sync_msk_1d(const int* lev, const int* n, const amrex_real* cf,
			   const int* plo, const int* phi, 
			   const int* lo, const int* hi);
  void mgt_set_sync_msk_2d(const int* lev, const int* n, const amrex_real* cf,
			   const int* plo, const int* phi, 
			   const int* lo, const int* hi);
  void mgt_set_sync_msk_3d(const int* lev, const int* n, const amrex_real* cf,
			   const int* plo, const int* phi, 
			   const int* lo, const int* hi);

  void mgt_set_vold_1d(const int* lev, const int* n, const amrex_real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  void mgt_set_vold_2d(const int* lev, const int* n, const amrex_real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  void mgt_set_vold_3d(const int* lev, const int* n, const amrex_real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);

  void mgt_get_sync_res_1d(const int* lev, const int* n, amrex_real* cf,
			   const int* plo, const int* phi, 
			   const int* lo, const int* hi);
  void mgt_get_sync_res_2d(const int* lev, const int* n, amrex_real* cf,
			   const int* plo, const int* phi, 
			   const int* lo, const int* hi);
  void mgt_get_sync_res_3d(const int* lev, const int* n, amrex_real* cf,
			   const int* plo, const int* phi, 
			   const int* lo, const int* hi);

  void mgt_set_rhcc_nodal_1d(const int* lev, const int* n, const amrex_real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  void mgt_set_rhcc_nodal_2d(const int* lev, const int* n, const amrex_real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);
  void mgt_set_rhcc_nodal_3d(const int* lev, const int* n, const amrex_real* cf,
		       const int* plo, const int* phi, 
		       const int* lo, const int* hi);

  void mgt_alloc_nodal_sync();
  void mgt_dealloc_nodal_sync();
  void mgt_compute_sync_resid_crse();
  void mgt_compute_sync_resid_fine();

  void mgt_alloc_rhcc_nodal();
  void mgt_dealloc_rhcc_nodal();
  void mgt_add_divucc();

  void mgt_dealloc();

  void mgt_nodal_dealloc();
  
  void mgt_solve(const amrex_real& tol, const amrex_real& abs_tol, const int* need_grad_phi, amrex_real* final_resnorm,
                 int* status, const int* always_use_bnorm);

  void mgt_applyop();
  
  void mgt_compute_flux(const int& n);
  
  void mgt_delete_flux(const int& n);
  
  void mgt_compute_residual();
  
  void mgt_nodal_solve(const amrex_real& tol, const amrex_real& abs_tol);
  
  void mgt_divu(int* lo_inflow, int* hi_inflow);
  
  void mgt_newu();

  void mgt_set_defaults(const int* nu1, const int* nu2, const int* nub, const int* nuf,
                        const int* max_iter, const int* bottom_max_iter,
                        const int* bottom_solver, const amrex_real* bottom_solver_eps,
                        const amrex_real* max_L0_growth,
                        const int* verbose, const int* cg_verbose,
                        const int* max_nlevel, const int* min_width,
                        const int* cycle, const int* smoother, const int* stencil_type);

  void mgt_set_nodal_defaults(const int* nu1, const int* nu2, const int* nub, const int* nuf,
                              const int* max_iter, const int* bottom_max_iter,
                              const int* bottom_solver, const amrex_real* bottom_solver_eps,
                              const int* verbose, const int* cg_verbose,
                              const int* max_nlevel, const int* min_width,
                              const int* cycle, const int* smoother, const int* stencil_type);

  void mgt_get_defaults(const int* nu1, const int* nu2, const int* nub, const int* nuf,
                        const int* max_iter, const int* bottom_max_iter,
                        const int* bottom_solver,
                        const amrex_real* max_L0_growth,
                        const int* verbose, const int* cg_verbose,
                        const int* max_nlevel, const int* min_width,
                        const int* cycle, const int* smoother);

  void mgt_get_nodal_defaults(const int* nu1, const int* nu2, const int* nub, const int* nuf,
                              const int* max_iter, const int* bottom_max_iter,
                              const int* bottom_solver,
                              const int* verbose, const int* cg_verbose,
                              const int* max_nlevel, const int* min_width,
                              const int* cycle, const int* smoother);

  void mgt_set_maxorder(const int* maxorder);

  int amrex_f90mg_get_niters ();
  
#ifdef __cplusplus
}
#endif

#endif
