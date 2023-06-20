#ifndef AMREX_MLLINOP_F_H_
#define AMREX_MLLINOP_F_H_
#include <AMReX_Config.H>

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_mllinop_apply_bc (const int* lo, const int* hi,
                                 amrex_real* phi, const int* philo, const int* phihi,
                                 const int* mask, const int* mlo, const int* mhi,
                                 int cdir, int bct, amrex_real bcl,
                                 const amrex_real* bcval, const int* blo, const int* bhi,
                                 int maxorder, const amrex_real* dxinv,
                                 int inhomog, int nc, int cross);

#ifdef __cplusplus
}
#endif

#endif
