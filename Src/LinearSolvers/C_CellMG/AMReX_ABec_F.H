#ifndef _AABBEC_F_H_
#define _AABBEC_F_H_

#include <AMReX_REAL.H>

#if        defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 1)
#define FORT_LINESOLVE     linesolve1daabbec
#define FORT_ADOTX         adotx1daabbec
#define FORT_NORMA         norma1daabbec
#define FORT_FLUX          flux1daabbec
#endif

#if (BL_SPACEDIM == 2)
#define FORT_GSRB          gsrb2daabbec
#define FORT_JACOBI        jacobi2daabbec
#define FORT_ADOTX         adotx2daabbec
#define FORT_NORMA         norma2daabbec
#define FORT_FLUX          flux2daabbec
#endif

#if (BL_SPACEDIM == 3)
#define FORT_GSRB          gsrb3daabbec
#define FORT_JACOBI        jacobi3daabbec
#define FORT_ADOTX         adotx3daabbec
#define FORT_NORMA         norma3daabbec
#define FORT_FLUX          flux3daabbec
#endif

#else

#if (BL_SPACEDIM == 1)
 
#if  defined(BL_FORT_USE_UPPERCASE)
#define FORT_LINESOLVE     LINESOLVE1DAABBEC
#define FORT_ADOTX    ADOTX1DAABBEC
#define FORT_NORMA    NORMA1DAABBEC
#define FORT_FLUX     FLUX1DAABBEC
#elif defined(BL_FORT_USE_LOWERCASE)
#define FORT_LINESOLVE     linesolve1daabbec_
#define FORT_ADOTX    adotx1daabbec
#define FORT_NORMA    norma1daabbec
#define FORT_FLUX     flux1daabbec
#elif defined(BL_FORT_USE_UNDERSCORE)
#define FORT_LINESOLVE     linesolve1daabbec_
#define FORT_ADOTX    adotx1daabbec_
#define FORT_NORMA    norma1daabbec_
#define FORT_FLUX     flux1daabbec_
#endif
 
#endif

#if (BL_SPACEDIM == 2)

#if  defined(BL_FORT_USE_UPPERCASE)
#define FORT_GSRB     GSRB2DAABBEC
#define FORT_JACOBI   JACOBI2DAABBEC
#define FORT_ADOTX    ADOTX2DAABBEC
#define FORT_NORMA    NORMA2DAABBEC
#define FORT_FLUX     FLUX2DAABBEC
#elif defined(BL_FORT_USE_LOWERCASE)
#define FORT_GSRB     gsrb2daabbec
#define FORT_JACOBI   jacobi2daabbec
#define FORT_ADOTX    adotx2daabbec
#define FORT_NORMA    norma2daabbec
#define FORT_FLUX     flux2daabbec
#elif defined(BL_FORT_USE_UNDERSCORE)
#define FORT_GSRB     gsrb2daabbec_
#define FORT_JACOBI   jacobi2daabbec_
#define FORT_ADOTX    adotx2daabbec_
#define FORT_NORMA    norma2daabbec_
#define FORT_FLUX     flux2daabbec_
#endif

#endif


#if (BL_SPACEDIM == 3)

#if   defined(BL_FORT_USE_UPPERCASE)
#define FORT_GSRB     GSRB3DAABBEC
#define FORT_JACOBI   JACOBI3DAABBEC
#define FORT_ADOTX    ADOTX3DAABBEC
#define FORT_NORMA    NORMA3DAABBEC
#define FORT_FLUX     FLUX3DAABBEC
#elif defined(BL_FORT_USE_LOWERCASE)
#define FORT_GSRB     gsrb3daabbec
#define FORT_JACOBI   jacobi3daabbec
#define FORT_ADOTX    adotx3daabbec
#define FORT_NORMA    norma3daabbec
#define FORT_FLUX     flux3daabbec
#elif defined(BL_FORT_USE_UNDERSCORE)
#define FORT_GSRB     gsrb3daabbec_
#define FORT_JACOBI   jacobi3daabbec_
#define FORT_ADOTX    adotx3daabbec_
#define FORT_NORMA    norma3daabbec_
#define FORT_FLUX     flux3daabbec_
#endif

#endif

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C"
{
#endif

#if (BL_SPACEDIM == 1)
    void FORT_LINESOLVE (
        amrex_real* phi       , ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const amrex_real* rhs , ARLIM_P(rhs_lo), ARLIM_P(phi_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a   , ARLIM_P(a_lo),   ARLIM_P(a_hi),
        const amrex_real* bX  , ARLIM_P(bX_lo),  ARLIM_P(bX_hi),
        const amrex_real* den0, ARLIM_P(den0_lo),ARLIM_P(den0_hi),
        const int* m0   , ARLIM_P(m0_lo),  ARLIM_P(m0_hi),
        const amrex_real* den2, ARLIM_P(den2_lo),ARLIM_P(den2_hi),
        const int* m2   , ARLIM_P(m2_lo),  ARLIM_P(m2_hi),
        const int* lo, const int* hi, const int *nc,
        const amrex_real *h
        );

    void FORT_ADOTX(
        amrex_real *y      , ARLIM_P(y_lo), ARLIM_P(y_hi),
        const amrex_real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const amrex_real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );
    
    void FORT_NORMA(
        amrex_real* res      ,
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const amrex_real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );
    
    void FORT_FLUX(
        const amrex_real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const amrex_real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const int *xlo, const int *xhi, const int *nc,
        const amrex_real *h,
        const amrex_real* xflux, ARLIM_P(xflux_lo), ARLIM_P(xflux_hi)
        );
#endif    
#if (BL_SPACEDIM == 2)
    void FORT_GSRB (
        amrex_real* phi       , ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const amrex_real* rhs , ARLIM_P(rhs_lo), ARLIM_P(phi_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a   , ARLIM_P(a_lo),   ARLIM_P(a_hi),
        const amrex_real* bX  , ARLIM_P(bX_lo),  ARLIM_P(bX_hi),
        const amrex_real* bY  , ARLIM_P(bY_lo),  ARLIM_P(bY_hi),
        const amrex_real* den0, ARLIM_P(den0_lo),ARLIM_P(den0_hi),
        const int* m0   , ARLIM_P(m0_lo),  ARLIM_P(m0_hi),
        const amrex_real* den1, ARLIM_P(den1_lo),ARLIM_P(den1_hi),
        const int* m1   , ARLIM_P(m1_lo),  ARLIM_P(m1_hi),
        const amrex_real* den2, ARLIM_P(den2_lo),ARLIM_P(den2_hi),
        const int* m2   , ARLIM_P(m2_lo),  ARLIM_P(m2_hi),
        const amrex_real* den3, ARLIM_P(den3_lo),ARLIM_P(den3_hi),
        const int* m3   , ARLIM_P(m3_lo),  ARLIM_P(m3_hi),
        const int* lo, const int* hi, const int* blo, const int* bhi, 
	const int *nc, const amrex_real *h, const  int* redblack
        );

    void FORT_JACOBI (
        amrex_real* phi       , ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const amrex_real* rhs , ARLIM_P(rhs_lo), ARLIM_P(phi_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a   , ARLIM_P(a_lo),   ARLIM_P(a_hi),
        const amrex_real* bX  , ARLIM_P(bX_lo),  ARLIM_P(bX_hi),
        const amrex_real* bY  , ARLIM_P(bY_lo),  ARLIM_P(bY_hi),
        const amrex_real* den0, ARLIM_P(den0_lo),ARLIM_P(den0_hi),
        const int* m0   , ARLIM_P(m0_lo),  ARLIM_P(m0_hi),
        const amrex_real* den1, ARLIM_P(den1_lo),ARLIM_P(den1_hi),
        const int* m1   , ARLIM_P(m1_lo),  ARLIM_P(m1_hi),
        const amrex_real* den2, ARLIM_P(den2_lo),ARLIM_P(den2_hi),
        const int* m2   , ARLIM_P(m2_lo),  ARLIM_P(m2_hi),
        const amrex_real* den3, ARLIM_P(den3_lo),ARLIM_P(den3_hi),
        const int* m3   , ARLIM_P(m3_lo),  ARLIM_P(m3_hi),
        const int* lo, const int* hi, const int *nc,
        const amrex_real *h
        );

    void FORT_ADOTX(
        amrex_real *y      , ARLIM_P(y_lo), ARLIM_P(y_hi),
        const amrex_real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const amrex_real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const amrex_real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );
    
    void FORT_NORMA(
        amrex_real* res      ,
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const amrex_real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const amrex_real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );
    
    void FORT_FLUX(
        const amrex_real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const amrex_real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const amrex_real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const int *xlo, const int *xhi, 
        const int *ylo, const int *yhi, 
	const int *nc,
        const amrex_real *h,
        const amrex_real* xflux, ARLIM_P(xflux_lo), ARLIM_P(xflux_hi),
        const amrex_real* yflux, ARLIM_P(yflux_lo), ARLIM_P(yflux_hi)
        );
#endif    

#if (BL_SPACEDIM == 3)
    void FORT_GSRB (
        amrex_real* phi,       ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const amrex_real* rhs, ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const amrex_real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const amrex_real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const amrex_real* bZ, ARLIM_P(bZ_lo), ARLIM_P(bZ_hi),
        const amrex_real* den0, ARLIM_P(den0_lo), ARLIM_P(den0_hi),
        const int* m0   , ARLIM_P(m0_lo),   ARLIM_P(m0_hi),
        const amrex_real* den1, ARLIM_P(den1_lo), ARLIM_P(den1_hi),
        const int* m1   , ARLIM_P(m1_lo),   ARLIM_P(m1_hi),
        const amrex_real* den2, ARLIM_P(den2_lo), ARLIM_P(den2_hi),
        const int* m2   , ARLIM_P(m2_lo),   ARLIM_P(m2_hi),
        const amrex_real* den3, ARLIM_P(den3_lo), ARLIM_P(den3_hi),
        const int* m3   , ARLIM_P(m3_lo),   ARLIM_P(m3_hi),
        const amrex_real* den4, ARLIM_P(den4_lo), ARLIM_P(den4_hi),
        const int* m4   , ARLIM_P(m4_lo),   ARLIM_P(m4_hi),
        const amrex_real* den5, ARLIM_P(den5_lo), ARLIM_P(den5_hi),
        const int* m5   , ARLIM_P(m5_lo),   ARLIM_P(m5_hi),
        const int* lo, const int* hi, const int* blo, const int* bhi, 
	const int *nc, const amrex_real *h, const  int* redblack
        );

    void FORT_JACOBI (
        amrex_real* phi,       ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const amrex_real* rhs, ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const amrex_real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const amrex_real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const amrex_real* bZ, ARLIM_P(bZ_lo), ARLIM_P(bZ_hi),
        const amrex_real* den0, ARLIM_P(den0_lo), ARLIM_P(den0_hi),
        const int* m0   , ARLIM_P(m0_lo),   ARLIM_P(m0_hi),
        const amrex_real* den1, ARLIM_P(den1_lo), ARLIM_P(den1_hi),
        const int* m1   , ARLIM_P(m1_lo),   ARLIM_P(m1_hi),
        const amrex_real* den2, ARLIM_P(den2_lo), ARLIM_P(den2_hi),
        const int* m2   , ARLIM_P(m2_lo),   ARLIM_P(m2_hi),
        const amrex_real* den3, ARLIM_P(den3_lo), ARLIM_P(den3_hi),
        const int* m3   , ARLIM_P(m3_lo),   ARLIM_P(m3_hi),
        const amrex_real* den4, ARLIM_P(den4_lo), ARLIM_P(den4_hi),
        const int* m4   , ARLIM_P(m4_lo),   ARLIM_P(m4_hi),
        const amrex_real* den5, ARLIM_P(den5_lo), ARLIM_P(den5_hi),
        const int* m5   , ARLIM_P(m5_lo),   ARLIM_P(m5_hi),
        const int* lo, const int* hi, const int *nc,
        const amrex_real *h
        );

    void FORT_ADOTX(
        amrex_real *y      , ARLIM_P(y_lo), ARLIM_P(y_hi),
        const amrex_real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const amrex_real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const amrex_real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const amrex_real* bZ, ARLIM_P(bZ_lo), ARLIM_P(bZ_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );
    
    void FORT_NORMA(
        amrex_real* res      ,
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const amrex_real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const amrex_real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const amrex_real* bZ, ARLIM_P(bZ_lo), ARLIM_P(bZ_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );
    
    void FORT_FLUX(
        const amrex_real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const amrex_real* alpha, const amrex_real* beta,
        const amrex_real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const amrex_real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const amrex_real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const amrex_real* bZ, ARLIM_P(bZ_lo), ARLIM_P(bZ_hi),
        const int *xlo, const int *xhi, 
        const int *ylo, const int *yhi, 
        const int *zlo, const int *zhi, 
	const int *nc,
        const amrex_real *h,
        amrex_real* xflux, ARLIM_P(xflux_lo), ARLIM_P(xflux_hi),
        amrex_real* yflux, ARLIM_P(yflux_lo), ARLIM_P(yflux_hi),
        amrex_real* zflux, ARLIM_P(zflux_lo), ARLIM_P(zflux_hi)
        );
#endif
#ifdef __cplusplus
}
#endif
#endif

#endif /*_AABBEC_F_H_*/
