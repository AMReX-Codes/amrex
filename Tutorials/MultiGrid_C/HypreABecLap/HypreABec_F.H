#ifndef _HYPREABEC_F_H_
#define _HYPREABEC_F_H_

#include <REAL.H>

#if        defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 2)
#define FORT_HPACOEF          hpacoef_2d
#define FORT_HPBCOEF          hpbcoef_2d
#define FORT_HPBMAT3          hpbmat3_2d
#define FORT_HPBVEC3          hpbvec3_2d
#elif (BL_SPACEDIM == 3)
#define FORT_HPACOEF          hpacoef_3d
#define FORT_HPBCOEF          hpbcoef_3d
#define FORT_HPBMAT3          hpbmat3_3d
#define FORT_HPBVEC3          hpbvec3_3d
#endif

#else

#if (BL_SPACEDIM == 2)
#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_HPACOEF          HPACOEF_2d 
#define FORT_HPBCOEF          HPBCOEF_2d 
#define FORT_HPBMAT3          HPBMAT3_2d
#define FORT_HPBVEC3          HPBVEC3_2d
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_HPACOEF          hpacoef_2d 
#define FORT_HPBCOEF          hpbcoef_2d 
#define FORT_HPBMAT3          hpbmat3_2d
#define FORT_HPBVEC3          hpbvec3_2d
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_HPACOEF          hpacoef_2d_
#define FORT_HPBCOEF          hpbcoef_2d_
#define FORT_HPBMAT3          hpbmat3_2d_
#define FORT_HPBVEC3          hpbvec3_2d_
#endif
#endif

#if (BL_SPACEDIM == 3)
#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_HPACOEF          HPACOEF_3d 
#define FORT_HPBCOEF          HPBCOEF_3d 
#define FORT_HPBMAT3          HPBMAT3_3d
#define FORT_HPBVEC3          HPBVEC3_3d
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_HPACOEF          hpacoef_3d 
#define FORT_HPBCOEF          hpbcoef_3d 
#define FORT_HPBMAT3          hpbmat3_3d
#define FORT_HPBVEC3          hpbvec3_3d
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_HPACOEF          hpacoef_3d_
#define FORT_HPBCOEF          hpbcoef_3d_ 
#define FORT_HPBMAT3          hpbmat3_3d_
#define FORT_HPBVEC3          hpbvec3_3d_
#endif
#endif

#include <ArrayLim.H>

extern "C" {
  void FORT_HPACOEF(Real* mat, Real* acoefs, ARLIM_P(alo), ARLIM_P(ahi),
		    const int* lo, const int* hi, const Real& sa);
  void FORT_HPBCOEF(Real* mat, Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
		    const int* lo, const int* hi, const Real& sb,
		    const Real* dx, const int& idim);
  void FORT_HPBMAT3(Real* mat, Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
		    const int* lo, const int* hi, const Real& sb, const Real* dx, 
		    const int& idim, const int& bct, const Real& bcl, 
		    const int* msk, ARLIM_P(mlo), ARLIM_P(mhi));
  void FORT_HPBVEC3(Real* vec, Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
		    const int* lo, const int* hi, const Real& sb, const Real* dx, 
		    const int& idim, const int& bct, const Real& bcl, 
		    const int* msk, ARLIM_P(mlo), ARLIM_P(mhi),
		    const Real* bcv, ARLIM_P(bvlo), ARLIM_P(bvhi));
};

#endif

#endif /* _HYPREABEC_F_H_ */
