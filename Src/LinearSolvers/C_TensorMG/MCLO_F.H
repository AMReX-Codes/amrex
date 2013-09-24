#ifndef _MCLO_F_H_
#define _MCLO_F_H_

#include <REAL.H>

#if defined(BL_LANG_FORT)
#  define FORT_AVERAGECC mcaveragecc
#  define FORT_AVERAGEEC mcaverageec
#  define FORT_HARMONIC_AVERAGEEC mcharaverageec
#  define FORT_RESIDL    mcresid
#  define FORT_APPLYBC   mcapplybc
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_AVERAGECC MCAVERAGECC
#    define FORT_AVERAGEEC MCAVERAGEEC
#    define FORT_HARMONIC_AVERAGEEC MCHARAVERAGEEC
#    define FORT_RESIDL    MCRESID
#    define FORT_APPLYBC   MCAPPLYBC
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_AVERAGECC mcaveragecc
#    define FORT_AVERAGEEC mcaverageec
#    define FORT_HARMONIC_AVERAGEEC mcharaverageec
#    define FORT_RESIDL    mcresid
#    define FORT_APPLYBC   mcapplybc
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_AVERAGECC mcaveragecc_
#    define FORT_AVERAGEEC mcaverageec_
#    define FORT_HARMONIC_AVERAGEEC mcharaverageec_
#    define FORT_RESIDL    mcresid_
#    define FORT_APPLYBC   mcapplybc_
#  endif

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
#if BL_SPACEDIM==2
		    const int *mask,   ARLIM_P(mask_lo),  ARLIM_P(mask_hi), 
		    const int *maskp,  ARLIM_P(maskp_lo),  ARLIM_P(maskp_hi), 
		    const int *maskm,  ARLIM_P(maskm_lo),  ARLIM_P(maskm_hi), 
#elif BL_SPACEDIM==3
		    const int *maskn,  ARLIM_P(maskn_lo),  ARLIM_P(maskn_hi), 
		    const int *maske,  ARLIM_P(maske_lo),  ARLIM_P(maske_hi), 
		    const int *maskw,  ARLIM_P(maskw_lo),  ARLIM_P(maskw_hi), 
		    const int *masks,  ARLIM_P(masks_lo),  ARLIM_P(masks_hi), 
		    const int *maskt,  ARLIM_P(maskt_lo),  ARLIM_P(maskt_hi), 
		    const int *maskb,  ARLIM_P(maskb_lo),  ARLIM_P(maskb_hi), 
#endif
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
#endif

#endif /*_MCLO_F_H_*/
