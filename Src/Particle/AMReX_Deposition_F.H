#ifndef AMREX_DEPOSITION_F_H__ 
#define AMREX_DEPOSITION_F_H__ 

#include <AMReX_BLFort.H>

extern "C"
{
  void deposit(const void *, int np, 
	       amrex::Real* rho, const int* lo, const int* hi,
	       const amrex::Real* plo, const amrex::Real* dx);
}

#endif
