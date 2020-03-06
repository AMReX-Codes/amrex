#ifndef TEST_REACT_F_H_
#define TEST_REACT_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void init_state(const int* lo, const int* hi,
                  amrex::Real* state, const int* s_lo, const int* s_hi,
                  const int* ncomp, const int* npts);
  void get_state(amrex::Real* state, const int* s_lo, const int* s_hi, const int* ncomp,
                 const int* i, const int* j, const int* k, const int* c,
                 amrex::Real* f);
  void set_state(amrex::Real* state, const int* s_lo, const int* s_hi, const int* ncomp,
                 const int* i, const int* j, const int* k, const int* c,
                 const amrex::Real* f);

#ifdef CUDA_CVODE
__device__  void cv_f_rhs_device(amrex::Real* y, amrex::Real* ydot);

__device__  void cv_f_jtv_device(amrex::Real* v, amrex::Real* Jv, amrex::Real* y, amrex::Real* ydot);
#endif

#ifdef __cplusplus
}
#endif

#endif
