#ifndef _AABBEC_F_H_
#define _AABBEC_F_H_

/*
** $Id: ABec_F.H,v 1.1 1998-03-24 07:05:51 almgren Exp $
*/

#include <REAL.H>

#if        defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 2)
#define FORT_GSRB          gsrb2daabbec
#define FORT_CGPRECND      cgprecnd2daabbec
#define FORT_ADOTX         adotx2daabbec
#endif

#if (BL_SPACEDIM == 3)
#define FORT_GSRB          gsrb3daabbec
#define FORT_CGPRECND      cgprecnd3daabbec
#define FORT_ADOTX         adotx3daabbec
#endif

#endif


#if        defined(BL_LANG_CC)

#if (BL_SPACEDIM == 2)

#ifdef BL_FORT_USE_UPPERCASE
#define FORT_GSRB     GSRB2DAABBEC
#define FORT_CGPRECND CGPRECND2DAABBEC
#define FORT_ADOTX    ADOTX2DAABBEC
#endif

#ifdef BL_FORT_USE_UNDERSCORE
#define FORT_GSRB     gsrb2daabbec_
#define FORT_CGPRECND cgprecnd2daabbec_
#define FORT_ADOTX    adotx2daabbec_
#endif

#endif


#if (BL_SPACEDIM == 3)

#ifdef BL_FORT_USE_UPPERCASE
#define FORT_GSRB     GSRB3DAABBEC
#define FORT_CGPRECND CGPRECND3DAABBEC
#define FORT_ADOTX    ADOTX3DAABBEC
#endif

#ifdef BL_FORT_USE_UNDERSCORE
#define FORT_GSRB     gsrb3daabbec_
#define FORT_CGPRECND cgprecnd3daabbec_
#define FORT_ADOTX    adotx3daabbec_
#endif

#endif

#include <ArrayLim.H>

extern "C"
{
#if (BL_SPACEDIM == 2)
    void FORT_GSRB (
        Real* phi       , ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const Real* rhs , ARLIM_P(rhs_lo), ARLIM_P(phi_hi),
        const Real* alpha, const Real* beta,
        const Real* a   , ARLIM_P(a_lo),   ARLIM_P(a_hi),
        const Real* bX  , ARLIM_P(bX_lo),  ARLIM_P(bX_hi),
        const Real* bY  , ARLIM_P(bY_lo),  ARLIM_P(bY_hi),
        const Real* den0, ARLIM_P(den0_lo),ARLIM_P(den0_hi),
        const int* m0   , ARLIM_P(m0_lo),  ARLIM_P(m0_hi),
        const Real* den1, ARLIM_P(den1_lo),ARLIM_P(den1_hi),
        const int* m1   , ARLIM_P(m1_lo),  ARLIM_P(m1_hi),
        const Real* den2, ARLIM_P(den2_lo),ARLIM_P(den2_hi),
        const int* m2   , ARLIM_P(m2_lo),  ARLIM_P(m2_hi),
        const Real* den3, ARLIM_P(den3_lo),ARLIM_P(den3_hi),
        const int* m3   , ARLIM_P(m3_lo),  ARLIM_P(m3_hi),
        const int* lo, const int* hi, const int *nc,
        const Real *h, const  int* redblack
        );

    void FORT_CGPRECND(
        Real *tmpz      , ARLIM_P(tmpz_lo), ARLIM_P(tmpz_hi),
        Real *rho,
        const Real *tmpr, ARLIM_P(tmpr_lo), ARLIM_P(tmpr_hi),
        const Real* alpha, const Real* beta,
        const Real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const Real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const Real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const int *lo, const int *hi, const int *nc,
        const Real *h
        );

    void FORT_ADOTX(
        Real *y      , ARLIM_P(y_lo), ARLIM_P(y_hi),
        const Real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const Real* alpha, const Real* beta,
        const Real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const Real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const Real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const int *lo, const int *hi, const int *nc,
        const Real *h
        );
#endif    

#if (BL_SPACEDIM == 3)
    void FORT_GSRB (
        Real* phi,       ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const Real* rhs, ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
        const Real* alpha, const Real* beta,
        const Real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const Real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const Real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const Real* bZ, ARLIM_P(bZ_lo), ARLIM_P(bZ_hi),
        const Real* den0, ARLIM_P(den0_lo), ARLIM_P(den0_hi),
        const int* m0   , ARLIM_P(m0_lo),   ARLIM_P(m0_hi),
        const Real* den1, ARLIM_P(den1_lo), ARLIM_P(den1_hi),
        const int* m1   , ARLIM_P(m1_lo),   ARLIM_P(m1_hi),
        const Real* den2, ARLIM_P(den2_lo), ARLIM_P(den2_hi),
        const int* m2   , ARLIM_P(m2_lo),   ARLIM_P(m2_hi),
        const Real* den3, ARLIM_P(den3_lo), ARLIM_P(den3_hi),
        const int* m3   , ARLIM_P(m3_lo),   ARLIM_P(m3_hi),
        const Real* den4, ARLIM_P(den4_lo), ARLIM_P(den4_hi),
        const int* m4   , ARLIM_P(m4_lo),   ARLIM_P(m4_hi),
        const Real* den5, ARLIM_P(den5_lo), ARLIM_P(den5_hi),
        const int* m5   , ARLIM_P(m5_lo),   ARLIM_P(m5_hi),
        const int* lo, const int* hi, const int *nc,
        const Real *h, const  int* redblack
        );

    void FORT_CGPRECND(
        Real *tmpz      , ARLIM_P(tmpz_lo), ARLIM_P(tmpz_hi),
        Real *rho,
        const Real *tmpr, ARLIM_P(tmpr_lo), ARLIM_P(tmpr_hi),
        const Real* alpha, const Real* beta,
        const Real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const Real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const Real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const Real* bZ, ARLIM_P(bZ_lo), ARLIM_P(bZ_hi),
        const int *lo, const int *hi, const int *nc,
        const Real *h
        );

    void FORT_ADOTX(
        Real *y      , ARLIM_P(y_lo), ARLIM_P(y_hi),
        const Real *x, ARLIM_P(x_lo), ARLIM_P(x_hi),
        const Real* alpha, const Real* beta,
        const Real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const Real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const Real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const Real* bZ, ARLIM_P(bZ_lo), ARLIM_P(bZ_hi),
        const int *lo, const int *hi, const int *nc,
        const Real *h
        );
#endif
}
#endif

#endif /*_AABBEC_F_H_*/
