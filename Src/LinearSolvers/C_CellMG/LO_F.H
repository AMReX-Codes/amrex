#ifndef _LO_F_H_
#define _LO_F_H_

/*
** $Id: LO_F.H,v 1.3 2002-04-26 22:10:02 lijewski Exp $
*/

#include <REAL.H>

#if        defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 2)
#define FORT_AVERAGECC          averagecc2dgen
#define FORT_AVERAGEEC          averageec2dgen
#define FORT_HARMONIC_AVERAGEEC haraverageec2dgen
#define FORT_APPLYBC            applybc2dgen
#define FORT_RESIDL             resid2dgen
#endif

#if (BL_SPACEDIM == 3)
#define FORT_AVERAGECC          averagecc3dgen
#define FORT_AVERAGEEC          averageec3dgen
#define FORT_HARMONIC_AVERAGEEC haraverageec3dgen
#define FORT_APPLYBC            applybc3dgen
#define FORT_RESIDL             resid3dgen
#endif

#else

#if (BL_SPACEDIM == 2)
#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_AVERAGECC          AVERAGECC2DGEN
#define FORT_AVERAGEEC          AVERAGEEC2DGEN
#define FORT_HARMONIC_AVERAGEEC HARAVERAGEEC2DGEN
#define FORT_APPLYBC            APPLYBC2DGEN
#define FORT_RESIDL             RESID2DGEN
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_AVERAGECC          averagecc2dgen
#define FORT_AVERAGEEC          averageec2dgen
#define FORT_HARMONIC_AVERAGEEC haraverageec2dgen
#define FORT_APPLYBC            applybc2dgen
#define FORT_RESIDL             resid2dgen
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_AVERAGECC          averagecc2dgen_
#define FORT_AVERAGEEC          averageec2dgen_
#define FORT_HARMONIC_AVERAGEEC haraverageec2dgen_
#define FORT_APPLYBC            applybc2dgen_
#define FORT_RESIDL             resid2dgen_
#endif
#endif

#if (BL_SPACEDIM == 3)
#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_AVERAGECC          AVERAGECC3DGEN
#define FORT_AVERAGEEC          AVERAGEEC3DGEN
#define FORT_HARMONIC_AVERAGEEC HARAVERAGEEC3DGEN
#define FORT_APPLYBC            APPLYBC3DGEN
#define FORT_RESIDL             RESID3DGEN
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_AVERAGECC          averagecc3dgen
#define FORT_AVERAGEEC          averageec3dgen
#define FORT_HARMONIC_AVERAGEEC haraverageec3dgen
#define FORT_APPLYBC            applybc3dgen
#define FORT_RESIDL             resid3dgen
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_AVERAGECC          averagecc3dgen_
#define FORT_AVERAGEEC          averageec3dgen_
#define FORT_HARMONIC_AVERAGEEC haraverageec3dgen_
#define FORT_APPLYBC            applybc3dgen_
#define FORT_RESIDL             resid3dgen_
#endif
#endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_RESIDL (
        Real* res      , ARLIM_P(res_lo), ARLIM_P(res_hi),
        const Real* rhs, ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
        const Real* phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const int* lo, const int* hi, const int* nc
        );

    void FORT_APPLYBC(
        const int *flagden, const int *flagbc, const int *maxorder,
        Real *phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const int *cdr,
        const int *bct,
        const Real *bcl,
        const Real *bcval, ARLIM_P(bcval_lo), ARLIM_P(bcval_hi),
        const int *mask,   ARLIM_P(mask_lo),  ARLIM_P(mask_hi), 
        Real *den,         ARLIM_P(den_lo),   ARLIM_P(den_hi),
        const int *lo, const int *hi, const int *nc,
        const Real *h
        );

    void FORT_AVERAGECC (
        Real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
        const Real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
        const int *tlo, const int *thi, const int *nc
        );

    void FORT_AVERAGEEC (
        Real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
        const Real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
        const int *tlo, const int *thi, const int *nc,
        const int *axis
        );

    void FORT_HARMONIC_AVERAGEEC (
        Real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
        const Real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
        const int *tlo, const int *thi, const int *nc,
        const int *axis
        );
}
#endif

#endif /*_LO_F_H_*/
