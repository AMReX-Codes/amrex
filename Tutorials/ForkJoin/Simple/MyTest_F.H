#ifndef MY_TEST_F_H_
#define MY_TEST_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void actual_init_poisson(const int* lo, const int* hi,
                             amrex_real* rhs, const int* rlo, const int* rhi,
                             amrex_real* exact, const int* elo, const int* ehi,
                             const amrex_real* problo, const amrex_real* probhi, const amrex_real* dx);

    void actual_init_abeclap(const int* lo, const int* hi, const int* glo, const int* ghi,
                             amrex_real* rhs, const int* rlo, const int* rhi,
                             amrex_real* exact, const int* elo, const int* ehi,
                             amrex_real* alpha, const int* alo, const int *ahi,
                             amrex_real* beta, const int* blo, const int* bhi,
                             amrex_real a, amrex_real b,
                             const amrex_real* problo, const amrex_real* probhi, const amrex_real* dx);

#ifdef __cplusplus
}
#endif

#endif
