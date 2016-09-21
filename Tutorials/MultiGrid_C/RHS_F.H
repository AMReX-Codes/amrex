#include <REAL.H>
#include <ArrayLim.H>

extern "C"
{
    void set_rhs(Real* rhs      , ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
		 const int* lo, const int* hi, const Real* dx, 
		 const Real& a, const Real& b, const Real& sigma, const Real& w,
		 const int& ibnd);
}

