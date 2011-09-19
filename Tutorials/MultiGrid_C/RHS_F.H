#ifndef _RHS_H_
#define _RHS_H_

#include <REAL.H>

#if        defined(BL_LANG_FORT)
#define FORT_SET_RHS            set_rhs
#else

#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_SET_RHS            SET_RHS
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_SET_RHS            set_rhs
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_SET_RHS            set_rhs_
#endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_SET_RHS (
        Real* rhs      , ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
        const int* lo, const int* hi, const Real* dx, 
	const Real& a, const Real& b, const Real& sigma, const Real& w,
	const int& ibnd);
}
#endif

#endif /*_RHS_F_H_*/
