#ifndef _MCLO_F_H_
#define _MCLO_F_H_

#include <REAL.H>

#if defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 2)
#define FORT_AVERAGECC mcaveragecc2dgen
#define FORT_AVERAGEEC mcaverageec2dgen
#define FORT_HARMONIC_AVERAGEEC mcharaverageec2dgen
#define FORT_RESIDL    mcresid2dgen
#endif

#if (BL_SPACEDIM == 3)
#define FORT_AVERAGECC mcaveragecc3dgen
#define FORT_AVERAGEEC mcaverageec3dgen
#define FORT_HARMONIC_AVERAGEEC mcharaverageec3dgen
#define FORT_RESIDL    mcresid3dgen
#endif

#endif /*defined(BL_LANG_FORT)*/


#if defined(BL_LANG_CC)

#if (BL_SPACEDIM == 2)
#ifdef BL_FORT_USE_UPPERCASE
#define FORT_AVERAGECC MCAVERAGECC2DGEN
#define FORT_AVERAGEEC MCAVERAGEEC2DGEN
#define FORT_HARMONIC_AVERAGEEC MCHARAVERAGEEC2DGEN
#define FORT_RESIDL    MCRESID2DGEN
#endif
#ifdef BL_FORT_USE_UNDERSCORE
#define FORT_AVERAGECC mcaveragecc2dgen_
#define FORT_AVERAGEEC mcaverageec2dgen_
#define FORT_HARMONIC_AVERAGEEC mcharaverageec2dgen_
#define FORT_RESIDL    mcresid2dgen_
#endif
#endif

#if (BL_SPACEDIM == 3)
#ifdef BL_FORT_USE_UPPERCASE
#define FORT_AVERAGECC MCAVERAGECC3DGEN
#define FORT_AVERAGEEC MCAVERAGEEC3DGEN
#define FORT_HARMONIC_AVERAGEEC MCHARAVERAGEEC3DGEN
#define FORT_RESIDL    MCRESID3DGEN
#endif
#ifdef BL_FORT_USE_UNDERSCORE
#define FORT_AVERAGECC mcaveragecc3dgen_
#define FORT_AVERAGEEC mcaverageec3dgen_
#define FORT_HARMONIC_AVERAGEEC mcharaverageec3dgen_
#define FORT_RESIDL    mcresid3dgen_
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
		    const int *flagden,  // 1 if want values in den
		    const int *flagbc,   // 1 for inhomogeneous
		    const int *maxorder, // maxorder of interpolant
		    Real *phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
		    const int *cdir,
		    const int *bct,
		    const Real *bcl,
		    const Real *bcval, ARLIM_P(bcval_lo), ARLIM_P(bcval_hi),
		    const int *mask,   ARLIM_P(mask_lo),  ARLIM_P(mask_hi), 
		    const int *maskp,  ARLIM_P(maskp_lo),  ARLIM_P(maskp_hi), 
		    const int *maskm,  ARLIM_P(maskm_lo),  ARLIM_P(maskm_hi), 
		    Real *den,         ARLIM_P(den_lo),   ARLIM_P(den_hi),
		    const Real *exttd, ARLIM_P(exttd_lo), ARLIM_P(exttd_hi),
		    Real *tander,      ARLIM_P(tander_lo),ARLIM_P(tander_hi),
		    const int *lo, const int *hi, 
		    const int *nc,
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
#endif /*defined(BL_LANG_CC)*/

#endif /*_MCLO_F_H_*/
