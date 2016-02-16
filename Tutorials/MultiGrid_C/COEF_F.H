#ifndef _COEFF_H_
#define _COEFF_H_

#include <REAL.H>

#if        defined(BL_LANG_FORT)

#define FORT_COMP_ASOL          comp_asol
#define FORT_SET_ALPHA          set_alpha
#define FORT_SET_CC_COEF        set_cc_coef

#else

#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_COMP_ASOL          COMP_ASOL
#define FORT_SET_ALPHA          SET_ALPHA
#define FORT_SET_CC_COEF        SET_CC_COEF
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_COMP_ASOL          comp_asol
#define FORT_SET_ALPHA          set_alpha
#define FORT_SET_CC_COEF        set_cc_coef
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_COMP_ASOL          comp_asol_
#define FORT_SET_ALPHA          set_alpha_
#define FORT_SET_CC_COEF        set_cc_coef_
#endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_COMP_ASOL (
        Real* asol      , ARLIM_P(asol_lo), ARLIM_P(asol_hi),
        const int* lo, const int* hi, const Real* dx, const int& ibnd,
        const Real* offset);

    void FORT_SET_ALPHA (
        Real* alpha     , ARLIM_P(alpha_lo), ARLIM_P(alpha_hi),
        const int* lo, const int* hi, const Real* dx);

    void FORT_SET_CC_COEF (
        Real* coef     , ARLIM_P(coef_lo), ARLIM_P(coef_hi),
        const int* lo, const int* hi, const Real* dx, const Real& sigma, const Real& w);
}
#endif

#endif /*_COEF_F_H_*/
