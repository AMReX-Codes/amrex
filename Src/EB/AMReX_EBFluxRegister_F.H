#ifndef AMREX_EBFLUXREGISTER_F_H_
#define AMREX_EBFLUXREGISTER_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_eb_flux_reg_crseadd_va (const int* lo, const int* hi,
                                       amrex_real* d, const int* dlo, const int* dhi,
                                       const int* flag, const int* fglo, const int* fghi,
                                       const amrex_real* fx, const int* fxlo, const int* fxhi,
#if (AMREX_SPACEDIM >= 2)
                                       const amrex_real* fy, const int* fylo, const int* fyhi,
#endif
#if (AMREX_SPACEDIM == 3)
                                       const amrex_real* fz, const int* fzlo, const int* fzhi,
#endif
                                       const amrex_real* vfrac, const int* vlo, const int* vhi,
                                       const amrex_real* ax, const int* axlo, const int* axhi,
#if (AMREX_SPACEDIM >= 2)
                                       const amrex_real* ay, const int* aylo, const int* ayhi,
#endif
#if (AMREX_SPACEDIM == 3)
                                       const amrex_real* az, const int* azlo, const int* azhi,
#endif
                                       const amrex_real* dx, const amrex_real* dt, const int* nc);


    void amrex_eb_flux_reg_fineadd_va (const int* lo, const int* hi,
                                       amrex_real* d, const int* dlo, const int* dhi,
                                       const amrex_real* f, const int* flo, const int* fhi,
                                       const amrex_real* cvol, const int* clo, const int* chi,
                                       const amrex_real* vfrac, const int* vlo, const int* vhi,
                                       const amrex_real* ax, const int* axlo, const int* axhi,
#if (AMREX_SPACEDIM >= 2)
                                       const amrex_real* ay, const int* aylo, const int* ayhi,
#endif
#if (AMREX_SPACEDIM == 3)
                                       const amrex_real* az, const int* azlo, const int* azhi,
#endif
                                       const amrex_real* dx, const amrex_real*dt, const int* nc,
                                       const int* dir, const int* side, const int* ratio);


    void amrex_eb_flux_reg_fineadd_dm(const int* lo, const int* hi,
                                      amrex_real* d, const int* dlo, const int* dhi,
                                      const amrex_real* dm, const int* mlo, const int* mhi,
                                      const amrex_real* cvol, const int* clo, const int* chi,
                                      const amrex_real* vfrac, const int* vlo, const int* vhi,
                                      const amrex_real* dx, const int* nc, const int* ratio);


    void amrex_eb_rereflux_from_crse (const int* lo, const int* hi,
                                      amrex_real* d, const int* dlo, const int* dhi,
                                      const amrex_real* s, const int* slo, const int* shi,
                                      const int* amrflg, const int* aflo, const int* afhi,
                                      const void* ebflg, const int* eflo, const int* efhi,
                                      const amrex_real* v, const int* vlo, const int* vhi,
                                      const int* nc);

    void amrex_eb_rereflux_to_fine (const int* lo, const int* hi,
                                    amrex_real* fd, const int* fdlo, const int* fdhi,
                                    const amrex_real* cd, const int* cdlo, const int* cdhi,
                                    const int* msk, const int * mlo, const int* mhi,
                                    const int* nc, const int* ratio);


#ifdef __cplusplus
}
#endif

#endif
