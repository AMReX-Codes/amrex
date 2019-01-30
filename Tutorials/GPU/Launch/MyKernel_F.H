#ifndef MY_KERNEL_F_H
#define MY_KERNEL_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

AMREX_GPU_DEVICE void plusone_cudafort (const int* lo, const int* hi,
                                        amrex_real* dat, const int* dlo, const int* dhi);

#ifdef AMREX_USE_ACC
void plusone_acc (const int* lo, const int* hi,
                  amrex_real* dat, const int* dlo, const int* dhi);
#endif

#ifdef AMREX_OMP_OFFLOAD
void plusone_omp (const int* lo, const int* hi,
                  amrex_real* dat, const int* dlo, const int* dhi);
#endif

#ifdef __cplusplus
}
#endif

#endif
