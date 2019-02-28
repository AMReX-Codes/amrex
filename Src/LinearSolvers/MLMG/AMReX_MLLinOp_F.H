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
                                 const int inhomog, const int nc, const int cross);

    void amrex_mllinop_comp_interp_coef0 (const int* lo, const int* hi,
                                          amrex_real* den, const int* dlo, const int* dhi,
                                          const int* mask, const int* mlo, const int* mhi,
                                          const int cdir, const int bct, const amrex_real bcl,
                                          const int maxorder, const amrex_real* dxinv, const int nc);

    void amrex_mllinop_apply_metric (const int* lo, const int* hi,
                                     amrex_real* d, const int* dlo, const int* dhi,
                                     const amrex_real* r, const int* rlo, const int* rhi,
				     const int nc);

    void amrex_mllinop_grad (const int* xlo, const int* xhi,
#if (AMREX_SPACEDIM >= 2)
                             const int* ylo, const int* yhi,
#if (AMREX_SPACEDIM == 3)
                             const int* zlo, const int* zhi,
#endif
#endif
                             const amrex_real* dat, const int* dlo, const int* dhi,
                             amrex_real* gx, const int* gxlo, const int* gxhi,
#if (AMREX_SPACEDIM >= 2)
                             amrex_real* gy, const int* gylo, const int* gyhi,
#if (AMREX_SPACEDIM == 3)
                             amrex_real* gz, const int* gzlo, const int* gzhi,
#endif
#endif
                             const amrex_real* dxinv);

#ifdef __cplusplus
}
#endif

#endif
