#ifndef MyTest_F_H_
#define MyTest_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif 

    void build_rhs(const int* lo, const int* hi, 
                   amrex_real* rhs, const int* rlo, const int* rhi, 
                   const amrex_real* a, const int* alo, const int* ahi, 
                   const amrex_real* bx, const int* bxlo, const int* bxhi, 
                   const amrex_real* by, const int* bylo, const int* byhi, 
                   const amrex_real* dx, const amrex_real* prob_lo, 
                   const amrex_real* prob_hi);

#ifdef __cplusplus
}
#endif
#endif 
