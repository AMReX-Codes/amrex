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
                               const amrex_real* sol, const int* slo, const int* shi,
#endif
                               const amrex_real* dxinv, const int face_only);


#ifdef __cplusplus
}
#endif

#endif

