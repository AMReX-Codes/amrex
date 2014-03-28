#ifndef _COEFF_H_
#define _COEFF_H_

#include <REAL.H>

#if        defined(BL_LANG_FORT)

#define FORT_COMP_ASOL          comp_asol
#define FORT_SET_ALPHA          set_alpha
#define FORT_SET_CC_COEF        set_cc_coef
#define FORT_COEF_TO_EDGES      coef_to_edges

#else

#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_COMP_ASOL          COMP_ASOL
#define FORT_SET_ALPHA          SET_ALPHA
#define FORT_SET_CC_COEF        SET_CC_COEF
#define FORT_COEF_TO_EDGES      COEF_TO_EDGES
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_COMP_ASOL          comp_asol
#define FORT_SET_ALPHA          set_alpha
#define FORT_SET_CC_COEF        set_cc_coef
#define FORT_COEF_TO_EDGES      coef_to_edges
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_COMP_ASOL          comp_asol_
#define FORT_SET_ALPHA          set_alpha_
#define FORT_SET_CC_COEF        set_cc_coef_
#define FORT_COEF_TO_EDGES      coef_to_edges_
#endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_COMP_ASOL (
        Real* asol      , ARLIM_P(asol_lo), ARLIM_P(asol_hi),
        const int* lo, const int* hi, const Real* dx, const int& ibnd);

    void FORT_SET_ALPHA (
        Real* alpha     , ARLIM_P(alpha_lo), ARLIM_P(alpha_hi),
        const int* lo, const int* hi, const Real* dx);

    void FORT_SET_CC_COEF (
        Real* coef     , ARLIM_P(coef_lo), ARLIM_P(coef_hi),
        const int* lo, const int* hi, const Real* dx, const Real& sigma, const Real& w);

    void FORT_COEF_TO_EDGES (int * dir, 
        Real* edge     , ARLIM_P(edge_lo), ARLIM_P(edge_hi),
        Real* coef     , ARLIM_P(coef_lo), ARLIM_P(coef_hi),
        const int* lo, const int* hi);
}
#endif

#endif /*_COEF_F_H_*/
