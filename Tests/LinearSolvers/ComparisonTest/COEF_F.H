#ifndef _COEFF_H_
#define _COEFF_H_

#include <AMReX_REAL.H>

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

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C"
{
#endif
  void FORT_SET_COEF(amrex_real* exac, amrex_real* alph, amrex_real* beta, amrex_real* rhs,
		     const int* lo, const int* hi, const amrex_real* prob_lo, const amrex_real* prob_hi, 
		     const amrex_real* dx, const amrex_real& a, const amrex_real& b, 
		     const amrex_real& sigma, const amrex_real& w, int& ibnd);
#ifdef __cplusplus
}
#endif
#endif

#endif /*_COEF_F_H_*/
