#ifndef _COEFF_H_
#define _COEFF_H_

#include <REAL.H>

#if        defined(BL_LANG_FORT)
#define FORT_SET_COEF        set_coef
#else

#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_SET_COEF        SET_COEF
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_SET_COEF        set_coef
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_SET_COEF        set_coef_
#endif

#include <ArrayLim.H>

extern "C"
{
  void FORT_SET_COEF(Real* exac, Real* alph, Real* beta, Real* rhs,
		     const int* lo, const int* hi, const Real* prob_lo, const Real* prob_hi, 
		     const Real* dx, const Real& a, const Real& b, 
		     const Real& sigma, const Real& w, int& ibnd);
}
#endif

#endif /*_COEF_F_H_*/
