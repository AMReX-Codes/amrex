#ifndef AMREX_MLPOISSON_F_H_
#define AMREX_MLPOISSON_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_mlpoisson_adotx (const int* lo, const int* hi,
                                amrex_real* y, const int* ylo, const int* yhi,
                                const amrex_real* x, const int* xlo, const int* xhi,
                                const amrex_real* dxinv);

    void amrex_mlpoisson_flux (const int* lo, const int* hi,
                               amrex_real* fx, const int* fxlo, const int* fxhi,
#if (AMREX_SPACEDIM >= 2)
                               amrex_real* fy, const int* fylo, const int* fyhi,
#if (AMREX_SPACEDIM == 3)
                               amrex_real* fz, const int* fzlo, const int* fzhi,
#endif
#endif
                               const amrex_real* sol, const int* slo, const int* shi,
                               const amrex_real* dxinv, const int face_only);

    void amrex_mlpoisson_gsrb (const int* lo, const int* hi,
                               amrex_real* phi, const int* philo, const int* phihi,
                               const amrex_real* rhs, const int* rlo, const int* rhi,
                               const amrex_real* f0, const int* f0lo, const int* f0hi,
                               const amrex_real* f1, const int* f1lo, const int* f1hi,
#if (AMREX_SPACEDIM >= 2)
                               const amrex_real* f2, const int* f2lo, const int* f2hi,
                               const amrex_real* f3, const int* f3lo, const int* f3hi,
#if (AMREX_SPACEDIM == 3)
                               const amrex_real* f4, const int* f4lo, const int* f4hi,
                               const amrex_real* f5, const int* f5lo, const int* f5hi,
#endif
#endif
                               const int* m0, const int* m0lo, const int* m0hi,
                               const int* m1, const int* m1lo, const int* m1hi,
#if (AMREX_SPACEDIM >= 2)
                               const int* m2, const int* m2lo, const int* m2hi,
                               const int* m3, const int* m3lo, const int* m3hi,
#if (AMREX_SPACEDIM == 3)
                               const int* m4, const int* m4lo, const int* m4hi,
                               const int* m5, const int* m5lo, const int* m5hi,
#endif
#endif
                               const int* blo, const int* bhi, const amrex_real* dxinv,
                               const int redblack);

#ifdef __cplusplus
}
#endif

#endif

