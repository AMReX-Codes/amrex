#ifndef AMREX_MLLINOP_F_H_
#define AMREX_MLLINOP_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_mllinop_apply_bc (const int* lo, const int* hi,
                                 amrex_real* phi, const int* philo, const int* phihi,
                                 const int* mask, const int* mlo, const int* mhi,
                                 const int cdir, const int bct, const amrex_real bcl,
                                 const amrex_real* bcval, const int* blo, const int* bhi,
                                 const int maxorder, const amrex_real* dxinv,
                                 const int inhomog);


    void amrex_mllinop_comp_interp_coef0 (const int* lo, const int* hi,
                                          amrex_real* den, const int* dlo, const int* dhi,
                                          const int* mask, const int* mlo, const int* mhi,
                                          const int cdir, const int bct, const amrex_real bcl,
                                          const int maxorder, const amrex_real* dxinv);

    void amrex_mllinop_apply_metric (const int* lo, const int* hi,
                                     amrex_real* d, const int* dlo, const int* dhi,
                                     const amrex_real* r, const int* rlo, const int* rhi);

    
#ifdef __cplusplus
}
#endif

#endif
