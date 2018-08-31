#ifndef MY_TEST_F_H_
#define MY_TEST_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void mytest_set_phi_reg(const int* lo, const int* hi,
                            const int* xlo, const int* xhi,
                            const int* ylo, const int* yhi,
#if (AMREX_SPACEDIM == 3)
                            const int* zlo, const int* zhi,
#endif
                            amrex_real* phie, const int* elo, const int* ehi,
                            amrex_real* rhs, const int* rlo, const int* rhi,
                            amrex_real* bx, const int* bxlo, const int* bxhi,
                            amrex_real* by, const int* bylo, const int* byhi,
#if (AMREX_SPACEDIM == 3)
                            amrex_real* bz, const int* bzlo, const int* bzhi,
#endif
                            const amrex_real* dx, const int* prob_type);


    void mytest_set_phi_eb(const int* lo, const int* hi,
                           const int* xlo, const int* xhi,
                           const int* ylo, const int* yhi,
#if (AMREX_SPACEDIM == 3)
                           const int* zlo, const int* zhi,
#endif
                           amrex_real* phie, const int* elo, const int* ehi,
                           amrex_real* phib, const int* blo, const int* bhi,
                           amrex_real* rhs, const int* rlo, const int* rhi,
                           amrex_real* bx, const int* bxlo, const int* bxhi,
                           amrex_real* by, const int* bylo, const int* byhi,
#if (AMREX_SPACEDIM == 3)
                           amrex_real* bz, const int* bzlo, const int* bzhi,
#endif
                           amrex_real* bb, const int* bblo, const int* bbhi,
                           const void* flag, const int* flo, const int* fhi,
                           const amrex_real* cent, const int* tlo, const int* thi,
                           const amrex_real* bcent, const int* clo, const int* chi,
                           const amrex_real* dx, const int* prob_type);


#ifdef __cplusplus
}
#endif

#endif

