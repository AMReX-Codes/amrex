#ifndef AMREX_YAFLUXREGISTER_F_H_
#define AMREX_YAFLUXREGISTER_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_ya_flux_reg_crseadd (const int* lo, const int* hi,
                                    amrex_real* d, const int* dlo, const int* dhi,
                                    const int* flag, const int* fglo, const int* fghi,
                                    const amrex_real* fx, const int* fxlo, const int* fxhi,
#if (AMREX_SPACEDIM >= 2)
                                    const amrex_real* fy, const int* fylo, const int* fyhi,
#endif
#if (AMREX_SPACEDIM == 3)
                                    const amrex_real* fz, const int* fzlo, const int* fzhi,
#endif
                                    const amrex_real* dx, const amrex_real* dt, const int* nc);


    void amrex_ya_flux_reg_fineadd (const int* lo, const int* hi,
                                    amrex_real* d, const int* dlo, const int* dhi,
                                    const amrex_real* f, const int* flo, const int* fhi,
                                    const amrex_real* dx, const amrex_real*dt, const int* nc,
                                    const int* dir, const int* side, const int* ratio);


#ifdef __cplusplus
}
#endif

#endif
