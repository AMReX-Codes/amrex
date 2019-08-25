#ifndef AMREX_EB2_F_H_
#define AMREX_EB2_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_eb2_build_cellflag_from_ap (const int* lo, const int*hi,
                                           const amrex_real* apx, const int* axlo, const int* axhi,
                                           const amrex_real* apy, const int* aylo, const int* ayhi,
#if (AMREX_SPACEDIM == 3)
                                           const amrex_real* apz, const int* azlo, const int* azhi,
#endif
                                           void* flag, const int* flo, const int* fhi);
    

#ifdef __cplusplus
}
#endif

#endif
