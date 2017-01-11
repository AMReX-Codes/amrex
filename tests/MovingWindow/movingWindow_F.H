#ifndef MYFUNC_F_H__ 
#define MUFUNC_F_H__ 

#include <BLFort.H>

extern "C"
{
  void init_E(Real* data, const int* lo, const int* hi, const int* ng, 
	      const Real* dx, const Real* prob_lo, const Real* prob_hi);

  void shift_E(Real* data, const int* lo, const int* hi, const int* ng, 
	       const Real* dx, const Real* prob_lo, const Real* prob_hi);

}
#endif
