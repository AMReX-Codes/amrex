#ifndef AMREX_EB2_F_H_
#define AMREX_EB2_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

#if (AMREX_SPACEDIM == 2)
    void amrex_eb2_gfab_build_types (const int* lo, const int* hi,
                                     const amrex_real* s, const int* slo, const int* shi,
                                     void* cell, const int* clo, const int* chi,
                                     void* fx, const int* fxlo, const int* fxhi,
                                     void* fy, const int* fylo, const int* fyhi,
                                     void* ex, const int* exlo, const int* exhi,
                                     void* ey, const int* eylo, const int* eyhi);
#elif (AMREX_SPACEDIM == 3)
    void amrex_eb2_gfab_build_types (const int* lo, const int* hi,
                                     const amrex_real* s, const int* slo, const int* shi,
                                     void* cell, const int* clo, const int* chi,
                                     void* fx, const int* fxlo, const int* fxhi,
                                     void* fy, const int* fylo, const int* fyhi,
                                     void* fz, const int* fzlo, const int* fzhi,
                                     void* ex, const int* exlo, const int* exhi,
                                     void* ey, const int* eylo, const int* eyhi,
                                     void* ez, const int* ezlo, const int* ezhi);
#endif

#ifdef __cplusplus
}
#endif

#endif
