#ifndef MOVINGWINDOW_F_H__ 
#define MOVINGWINDOW_F_H__ 

#include <AMReX_BLFort.H>

extern "C"
{
  void init_E(amrex::Real* data, const int* lo, const int* hi, const int* ng, 
	      const amrex::Real* dx, const amrex::Real* prob_lo, const amrex::Real* prob_hi);

}
#endif
