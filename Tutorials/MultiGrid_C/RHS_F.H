#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void set_rhs(amrex_real* rhs      , ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
		 const int* lo, const int* hi, const amrex_real* dx, 
		 const amrex_real& a, const amrex_real& b, const amrex_real& sigma, const amrex_real& w,
		 const int& ibnd);
#ifdef __cplusplus
}
#endif

