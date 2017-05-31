#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void comp_asol(amrex_real* asol      , ARLIM_P(asol_lo), ARLIM_P(asol_hi),
		    const int* lo, const int* hi, const amrex_real* dx, const int& ibnd,
		    const amrex_real* offset);

    void set_alpha(amrex_real* alpha     , ARLIM_P(alpha_lo), ARLIM_P(alpha_hi),
		   const int* lo, const int* hi, const amrex_real* dx);

    void set_cc_coef(amrex_real* coef     , ARLIM_P(coef_lo), ARLIM_P(coef_hi),
		     const int* lo, const int* hi, const amrex_real* dx, const amrex_real& sigma, const amrex_real& w);
#ifdef __cplusplus
}
#endif


