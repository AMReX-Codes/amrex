#ifndef _HYPREABEC_F_H_
#define _HYPREABEC_F_H_

#include <REAL.H>

#if        defined(BL_LANG_FORT)

#define FORT_HPACOEF          hpacoef
#define FORT_HPBCOEF          hpbcoef
#define FORT_HPMAT            hpmat
#define FORT_HPMAT3           hpmat3
#define FORT_HPBVEC           hpbvec
#define FORT_HPBVEC3          hpbvec3

#else

#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_HPACOEF          HPACOEF
#define FORT_HPBCOEF          HPBCOEF
#define FORT_HPMAT            HPMAT
#define FORT_HPMAT3           HPMAT3
#define FORT_HPBVEC           HPBVEC
#define FORT_HPBVEC3          HPBVEC3
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_HPACOEF          hpacoef
#define FORT_HPBCOEF          hpbcoef
#define FORT_HPMAT            hpmat
#define FORT_HPMAT3           hpmat3
#define FORT_HPBVEC           hpbvec
#define FORT_HPBVEC3          hpbvec3
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_HPACOEF          hpacoef_
#define FORT_HPBCOEF          hpbcoef_ 
#define FORT_HPMAT            hpmat_
#define FORT_HPMAT3           hpmat3_
#define FORT_HPBVEC           hpbvec_
#define FORT_HPBVEC3          hpbvec3_
#endif

#include <ArrayLim.H>

extern "C" {
  void FORT_HPACOEF(Real* mat, Real* acoefs, ARLIM_P(alo), ARLIM_P(ahi),
		    const int* lo, const int* hi, const Real& sa);
  void FORT_HPBCOEF(Real* mat, Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
		    const int* lo, const int* hi, const Real& sb,
		    const Real* dx, const int& idim);
  void FORT_HPMAT  (Real* mat, Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
		    const int* lo, const int* hi, const Real& sb, const Real* dx, 
		    const int& idim, const int& bct, const Real& bcl, 
		    const int* msk, ARLIM_P(mlo), ARLIM_P(mhi));
  void FORT_HPMAT3 (Real* mat, Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
		    const int* lo, const int* hi, const Real& sb, const Real* dx, 
		    const int& idim, const int& bct, const Real& bcl, 
		    const int* msk, ARLIM_P(mlo), ARLIM_P(mhi));
  void FORT_HPBVEC3(Real* vec, Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
		    const int* lo, const int* hi, const Real& sb, const Real* dx, 
		    const int& idim, const int& bct, const Real& bcl, 
		    const int* msk, ARLIM_P(mlo), ARLIM_P(mhi),
		    const Real* bcv, ARLIM_P(bvlo), ARLIM_P(bvhi));
  void FORT_HPBVEC (Real* vec, Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
		    const int* lo, const int* hi, const Real& sb, const Real* dx, 
		    const int& idim, const int& bct, const Real& bcl, 
		    const int* msk, ARLIM_P(mlo), ARLIM_P(mhi),
		    const Real* bcv, ARLIM_P(bvlo), ARLIM_P(bvhi));
};

#endif

#endif /* _HYPREABEC_F_H_ */
