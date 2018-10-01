#ifndef MY_TEST_F_H_
#define MY_TEST_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void init_prob (const int* lo, const int* hi,
                    amrex_real* rhs, const int* rlo, const int* rhi,
                    amrex_real* phi, const int* hlo, const int* hhi,
                    const amrex_real* dx);

#ifdef __cplusplus
}
#endif

#endif
