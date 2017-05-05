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
  void FORT_SET_COEF(amrex::Real* exac, amrex::Real* alph, amrex::Real* beta, amrex::Real* rhs,
		     const int* lo, const int* hi, const amrex::Real* prob_lo, const amrex::Real* prob_hi, 
		     const amrex::Real* dx, const amrex::Real& a, const amrex::Real& b, 
		     const amrex::Real& sigma, const amrex::Real& w, int& ibnd);
#ifdef __cplusplus
}
#endif
#endif

#endif /*_COEF_F_H_*/
