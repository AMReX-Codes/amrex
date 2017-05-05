#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void set_rhs(amrex::Real* rhs      , ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
		 const int* lo, const int* hi, const amrex::Real* dx, 
		 const amrex::Real& a, const amrex::Real& b, const amrex::Real& sigma, const amrex::Real& w,
		 const int& ibnd);
#ifdef __cplusplus
}
#endif

