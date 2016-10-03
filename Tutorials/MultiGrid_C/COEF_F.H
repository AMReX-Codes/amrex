#include <REAL.H>
#include <ArrayLim.H>

extern "C"
{
    void comp_asol(Real* asol      , ARLIM_P(asol_lo), ARLIM_P(asol_hi),
		    const int* lo, const int* hi, const Real* dx, const int& ibnd,
		    const Real* offset);

    void set_alpha(Real* alpha     , ARLIM_P(alpha_lo), ARLIM_P(alpha_hi),
		   const int* lo, const int* hi, const Real* dx);

    void set_cc_coef(Real* coef     , ARLIM_P(coef_lo), ARLIM_P(coef_hi),
		     const int* lo, const int* hi, const Real* dx, const Real& sigma, const Real& w);
}


