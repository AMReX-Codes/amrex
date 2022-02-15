#ifndef MY_FORT_H_
#define MY_FORT_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void fort_set_coef (const int* lo, const int* hi,
                        amrex_real* exact, const int* elo, const int* ehi,
                        amrex_real* alpha, const int* alo, const int* ahi,
                        amrex_real* beta,  const int* blo, const int* bhi,
                        amrex_real* rhs,   const int* rlo, const int* rhi,
                        const amrex_real* dx, const amrex_real* problo, const amrex_real* probhi,
                        const amrex_real* a, const amrex_real* b,
                        const amrex_real* sigma, const amrex_real* w, const char* bct);

#ifdef __cplusplus
}
#endif

#endif
