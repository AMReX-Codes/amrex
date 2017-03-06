#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>

extern "C"
{
    void comp_asol(amrex::Real* asol      , ARLIM_P(asol_lo), ARLIM_P(asol_hi),
		    const int* lo, const int* hi, const amrex::Real* dx, const int& ibnd,
		    const amrex::Real* offset);

    void set_alpha(amrex::Real* alpha     , ARLIM_P(alpha_lo), ARLIM_P(alpha_hi),
		   const int* lo, const int* hi, const amrex::Real* dx);

    void set_cc_coef(amrex::Real* coef     , ARLIM_P(coef_lo), ARLIM_P(coef_hi),
		     const int* lo, const int* hi, const amrex::Real* dx, const amrex::Real& sigma, const amrex::Real& w);
}


