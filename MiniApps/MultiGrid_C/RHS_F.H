#ifndef _RHS_H_
#define _RHS_H_

#include <AMReX_REAL.H>

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

#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void FORT_SET_RHS (
        amrex_real* rhs      , ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
        const int* lo, const int* hi, const amrex_real* dx, 
	const amrex_real& a, const amrex_real& b, const amrex_real& sigma, const amrex_real& w,
	const int& ibnd);
#ifdef __cplusplus
}
#endif
#endif

#endif /*_RHS_F_H_*/
