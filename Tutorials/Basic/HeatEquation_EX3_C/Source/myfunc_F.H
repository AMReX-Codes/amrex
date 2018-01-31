#ifndef MYFUNC_F_H__ 
#define MUFUNC_F_H__ 

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void init_phi(const int* lo, const int* hi,
                  amrex_real* data, const int* dlo, const int* dhi,
                  const amrex_real* dx, const amrex_real* prob_lo, const amrex_real* prob_hi);
    
#ifdef __cplusplus
}
#endif

#endif
