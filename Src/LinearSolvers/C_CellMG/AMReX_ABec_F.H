#ifndef AMREX_AABBEC_F_H_
#define AMREX_AABBEC_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

#if (BL_SPACEDIM == 1)
    void amrex_abec_linesolve (
        amrex_real* phi       , AMREX_ARLIM_P(phi_lo), AMREX_ARLIM_P(phi_hi),
        const amrex_real* rhs , AMREX_ARLIM_P(rhs_lo), AMREX_ARLIM_P(phi_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a   , AMREX_ARLIM_P(a_lo),   AMREX_ARLIM_P(a_hi),
        const amrex_real* bX  , AMREX_ARLIM_P(bX_lo),  AMREX_ARLIM_P(bX_hi),
        const amrex_real* den0, AMREX_ARLIM_P(den0_lo),AMREX_ARLIM_P(den0_hi),
        const int* m0   , AMREX_ARLIM_P(m0_lo),  AMREX_ARLIM_P(m0_hi),
        const amrex_real* den2, AMREX_ARLIM_P(den2_lo),AMREX_ARLIM_P(den2_hi),
        const int* m2   , AMREX_ARLIM_P(m2_lo),  AMREX_ARLIM_P(m2_hi),
        const int* lo, const int* hi, const int *nc,
        const amrex_real *h
        );

    void amrex_abec_adotx(
        amrex_real *y      , AMREX_ARLIM_P(y_lo), AMREX_ARLIM_P(y_hi),
        const amrex_real *x, AMREX_ARLIM_P(x_lo), AMREX_ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , AMREX_ARLIM_P(a_lo),  AMREX_ARLIM_P(a_hi),
        const amrex_real* bX, AMREX_ARLIM_P(bX_lo), AMREX_ARLIM_P(bX_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );
    
    void amrex_abec_norma(
        amrex_real* res      ,
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , AMREX_ARLIM_P(a_lo),  AMREX_ARLIM_P(a_hi),
        const amrex_real* bX, AMREX_ARLIM_P(bX_lo), AMREX_ARLIM_P(bX_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );
    
    void amrex_abec_flux(
        const amrex_real *x, AMREX_ARLIM_P(x_lo), AMREX_ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , AMREX_ARLIM_P(a_lo),  AMREX_ARLIM_P(a_hi),
        const amrex_real* bX, AMREX_ARLIM_P(bX_lo), AMREX_ARLIM_P(bX_hi),
        const int *xlo, const int *xhi, const int *nc,
        const amrex_real *h,
        const amrex_real* xflux, AMREX_ARLIM_P(xflux_lo), AMREX_ARLIM_P(xflux_hi)
        );
#endif    
#if (BL_SPACEDIM == 2)
    void amrex_abec_gsrb (
        amrex_real* phi       , AMREX_ARLIM_P(phi_lo), AMREX_ARLIM_P(phi_hi),
        const amrex_real* rhs , AMREX_ARLIM_P(rhs_lo), AMREX_ARLIM_P(phi_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a   , AMREX_ARLIM_P(a_lo),   AMREX_ARLIM_P(a_hi),
        const amrex_real* bX  , AMREX_ARLIM_P(bX_lo),  AMREX_ARLIM_P(bX_hi),
        const amrex_real* bY  , AMREX_ARLIM_P(bY_lo),  AMREX_ARLIM_P(bY_hi),
        const amrex_real* den0, AMREX_ARLIM_P(den0_lo),AMREX_ARLIM_P(den0_hi),
        const int* m0   , AMREX_ARLIM_P(m0_lo),  AMREX_ARLIM_P(m0_hi),
        const amrex_real* den1, AMREX_ARLIM_P(den1_lo),AMREX_ARLIM_P(den1_hi),
        const int* m1   , AMREX_ARLIM_P(m1_lo),  AMREX_ARLIM_P(m1_hi),
        const amrex_real* den2, AMREX_ARLIM_P(den2_lo),AMREX_ARLIM_P(den2_hi),
        const int* m2   , AMREX_ARLIM_P(m2_lo),  AMREX_ARLIM_P(m2_hi),
        const amrex_real* den3, AMREX_ARLIM_P(den3_lo),AMREX_ARLIM_P(den3_hi),
        const int* m3   , AMREX_ARLIM_P(m3_lo),  AMREX_ARLIM_P(m3_hi),
        const int* lo, const int* hi, const int* blo, const int* bhi, 
	const int *nc, const amrex_real *h, const  int* redblack
        );

    void amrex_abec_jacobi (
        amrex_real* phi       , AMREX_ARLIM_P(phi_lo), AMREX_ARLIM_P(phi_hi),
        const amrex_real* rhs , AMREX_ARLIM_P(rhs_lo), AMREX_ARLIM_P(phi_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a   , AMREX_ARLIM_P(a_lo),   AMREX_ARLIM_P(a_hi),
        const amrex_real* bX  , AMREX_ARLIM_P(bX_lo),  AMREX_ARLIM_P(bX_hi),
        const amrex_real* bY  , AMREX_ARLIM_P(bY_lo),  AMREX_ARLIM_P(bY_hi),
        const amrex_real* den0, AMREX_ARLIM_P(den0_lo),AMREX_ARLIM_P(den0_hi),
        const int* m0   , AMREX_ARLIM_P(m0_lo),  AMREX_ARLIM_P(m0_hi),
        const amrex_real* den1, AMREX_ARLIM_P(den1_lo),AMREX_ARLIM_P(den1_hi),
        const int* m1   , AMREX_ARLIM_P(m1_lo),  AMREX_ARLIM_P(m1_hi),
        const amrex_real* den2, AMREX_ARLIM_P(den2_lo),AMREX_ARLIM_P(den2_hi),
        const int* m2   , AMREX_ARLIM_P(m2_lo),  AMREX_ARLIM_P(m2_hi),
        const amrex_real* den3, AMREX_ARLIM_P(den3_lo),AMREX_ARLIM_P(den3_hi),
        const int* m3   , AMREX_ARLIM_P(m3_lo),  AMREX_ARLIM_P(m3_hi),
        const int* lo, const int* hi, const int *nc,
        const amrex_real *h
        );

    void amrex_abec_adotx(
        amrex_real *y      , AMREX_ARLIM_P(y_lo), AMREX_ARLIM_P(y_hi),
        const amrex_real *x, AMREX_ARLIM_P(x_lo), AMREX_ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , AMREX_ARLIM_P(a_lo),  AMREX_ARLIM_P(a_hi),
        const amrex_real* bX, AMREX_ARLIM_P(bX_lo), AMREX_ARLIM_P(bX_hi),
        const amrex_real* bY, AMREX_ARLIM_P(bY_lo), AMREX_ARLIM_P(bY_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );
    
    void amrex_abec_norma(
        amrex_real* res      ,
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , AMREX_ARLIM_P(a_lo),  AMREX_ARLIM_P(a_hi),
        const amrex_real* bX, AMREX_ARLIM_P(bX_lo), AMREX_ARLIM_P(bX_hi),
        const amrex_real* bY, AMREX_ARLIM_P(bY_lo), AMREX_ARLIM_P(bY_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );
    
    void amrex_abec_flux(
        const amrex_real *x, AMREX_ARLIM_P(x_lo), AMREX_ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , AMREX_ARLIM_P(a_lo),  AMREX_ARLIM_P(a_hi),
        const amrex_real* bX, AMREX_ARLIM_P(bX_lo), AMREX_ARLIM_P(bX_hi),
        const amrex_real* bY, AMREX_ARLIM_P(bY_lo), AMREX_ARLIM_P(bY_hi),
        const int *xlo, const int *xhi, 
        const int *ylo, const int *yhi, 
	const int *nc,
        const amrex_real *h,
        const amrex_real* xflux, AMREX_ARLIM_P(xflux_lo), AMREX_ARLIM_P(xflux_hi),
        const amrex_real* yflux, AMREX_ARLIM_P(yflux_lo), AMREX_ARLIM_P(yflux_hi)
        );
#endif    

#if (BL_SPACEDIM == 3)
    void amrex_abec_gsrb (
        amrex_real* phi,       AMREX_ARLIM_P(phi_lo), AMREX_ARLIM_P(phi_hi),
        const amrex_real* rhs, AMREX_ARLIM_P(rhs_lo), AMREX_ARLIM_P(rhs_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , AMREX_ARLIM_P(a_lo),  AMREX_ARLIM_P(a_hi),
        const amrex_real* bX, AMREX_ARLIM_P(bX_lo), AMREX_ARLIM_P(bX_hi),
        const amrex_real* bY, AMREX_ARLIM_P(bY_lo), AMREX_ARLIM_P(bY_hi),
        const amrex_real* bZ, AMREX_ARLIM_P(bZ_lo), AMREX_ARLIM_P(bZ_hi),
        const amrex_real* den0, AMREX_ARLIM_P(den0_lo), AMREX_ARLIM_P(den0_hi),
        const int* m0   , AMREX_ARLIM_P(m0_lo),   AMREX_ARLIM_P(m0_hi),
        const amrex_real* den1, AMREX_ARLIM_P(den1_lo), AMREX_ARLIM_P(den1_hi),
        const int* m1   , AMREX_ARLIM_P(m1_lo),   AMREX_ARLIM_P(m1_hi),
        const amrex_real* den2, AMREX_ARLIM_P(den2_lo), AMREX_ARLIM_P(den2_hi),
        const int* m2   , AMREX_ARLIM_P(m2_lo),   AMREX_ARLIM_P(m2_hi),
        const amrex_real* den3, AMREX_ARLIM_P(den3_lo), AMREX_ARLIM_P(den3_hi),
        const int* m3   , AMREX_ARLIM_P(m3_lo),   AMREX_ARLIM_P(m3_hi),
        const amrex_real* den4, AMREX_ARLIM_P(den4_lo), AMREX_ARLIM_P(den4_hi),
        const int* m4   , AMREX_ARLIM_P(m4_lo),   AMREX_ARLIM_P(m4_hi),
        const amrex_real* den5, AMREX_ARLIM_P(den5_lo), AMREX_ARLIM_P(den5_hi),
        const int* m5   , AMREX_ARLIM_P(m5_lo),   AMREX_ARLIM_P(m5_hi),
        const int* lo, const int* hi, const int* blo, const int* bhi, 
	const int *nc, const amrex_real *h, const  int* redblack
        );

    void amrex_abec_jacobi (
        amrex_real* phi,       AMREX_ARLIM_P(phi_lo), AMREX_ARLIM_P(phi_hi),
        const amrex_real* rhs, AMREX_ARLIM_P(rhs_lo), AMREX_ARLIM_P(rhs_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , AMREX_ARLIM_P(a_lo),  AMREX_ARLIM_P(a_hi),
        const amrex_real* bX, AMREX_ARLIM_P(bX_lo), AMREX_ARLIM_P(bX_hi),
        const amrex_real* bY, AMREX_ARLIM_P(bY_lo), AMREX_ARLIM_P(bY_hi),
        const amrex_real* bZ, AMREX_ARLIM_P(bZ_lo), AMREX_ARLIM_P(bZ_hi),
        const amrex_real* den0, AMREX_ARLIM_P(den0_lo), AMREX_ARLIM_P(den0_hi),
        const int* m0   , AMREX_ARLIM_P(m0_lo),   AMREX_ARLIM_P(m0_hi),
        const amrex_real* den1, AMREX_ARLIM_P(den1_lo), AMREX_ARLIM_P(den1_hi),
        const int* m1   , AMREX_ARLIM_P(m1_lo),   AMREX_ARLIM_P(m1_hi),
        const amrex_real* den2, AMREX_ARLIM_P(den2_lo), AMREX_ARLIM_P(den2_hi),
        const int* m2   , AMREX_ARLIM_P(m2_lo),   AMREX_ARLIM_P(m2_hi),
        const amrex_real* den3, AMREX_ARLIM_P(den3_lo), AMREX_ARLIM_P(den3_hi),
        const int* m3   , AMREX_ARLIM_P(m3_lo),   AMREX_ARLIM_P(m3_hi),
        const amrex_real* den4, AMREX_ARLIM_P(den4_lo), AMREX_ARLIM_P(den4_hi),
        const int* m4   , AMREX_ARLIM_P(m4_lo),   AMREX_ARLIM_P(m4_hi),
        const amrex_real* den5, AMREX_ARLIM_P(den5_lo), AMREX_ARLIM_P(den5_hi),
        const int* m5   , AMREX_ARLIM_P(m5_lo),   AMREX_ARLIM_P(m5_hi),
        const int* lo, const int* hi, const int *nc,
        const amrex_real *h
        );

    void amrex_abec_adotx(
        amrex_real *y      , AMREX_ARLIM_P(y_lo), AMREX_ARLIM_P(y_hi),
        const amrex_real *x, AMREX_ARLIM_P(x_lo), AMREX_ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , AMREX_ARLIM_P(a_lo),  AMREX_ARLIM_P(a_hi),
        const amrex_real* bX, AMREX_ARLIM_P(bX_lo), AMREX_ARLIM_P(bX_hi),
        const amrex_real* bY, AMREX_ARLIM_P(bY_lo), AMREX_ARLIM_P(bY_hi),
        const amrex_real* bZ, AMREX_ARLIM_P(bZ_lo), AMREX_ARLIM_P(bZ_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );
    
    void amrex_abec_norma(
        amrex_real* res      ,
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , AMREX_ARLIM_P(a_lo),  AMREX_ARLIM_P(a_hi),
        const amrex_real* bX, AMREX_ARLIM_P(bX_lo), AMREX_ARLIM_P(bX_hi),
        const amrex_real* bY, AMREX_ARLIM_P(bY_lo), AMREX_ARLIM_P(bY_hi),
        const amrex_real* bZ, AMREX_ARLIM_P(bZ_lo), AMREX_ARLIM_P(bZ_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );
    
    void amrex_abec_flux(
        const amrex_real *x, AMREX_ARLIM_P(x_lo), AMREX_ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , AMREX_ARLIM_P(a_lo),  AMREX_ARLIM_P(a_hi),
        const amrex_real* bX, AMREX_ARLIM_P(bX_lo), AMREX_ARLIM_P(bX_hi),
        const amrex_real* bY, AMREX_ARLIM_P(bY_lo), AMREX_ARLIM_P(bY_hi),
        const amrex_real* bZ, AMREX_ARLIM_P(bZ_lo), AMREX_ARLIM_P(bZ_hi),
        const int *xlo, const int *xhi, 
        const int *ylo, const int *yhi, 
        const int *zlo, const int *zhi, 
	const int *nc,
        const amrex_real *h,
        amrex_real* xflux, AMREX_ARLIM_P(xflux_lo), AMREX_ARLIM_P(xflux_hi),
        amrex_real* yflux, AMREX_ARLIM_P(yflux_lo), AMREX_ARLIM_P(yflux_hi),
        amrex_real* zflux, AMREX_ARLIM_P(zflux_lo), AMREX_ARLIM_P(zflux_hi)
        );
#endif
#ifdef __cplusplus
}
#endif


#endif /*_AABBEC_F_H_*/
