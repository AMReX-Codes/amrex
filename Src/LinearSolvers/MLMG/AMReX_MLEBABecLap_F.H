#ifndef AMREX_MLEBABECLAP_F_H_
#define AMREX_MLEBABECLAP_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_mlebabeclap_adotx (const int* lo, const int* hi,
                                  amrex_real* y, const int* ylo, const int* yhi,
                                  const amrex_real* x, const int* xlo, const int* xhi,
                                  const amrex_real* a, const int* alo, const int* ahi,
                                  const amrex_real* bx, const int* bxlo, const int* bxhi,
#if (AMREX_SPACEDIM >= 2)
                                  const amrex_real* by, const int* bylo, const int* byhi,
#if (AMREX_SPACEDIM == 3)
                                  const amrex_real* bz, const int* bzlo, const int* bzhi,
#endif
#endif
                                  const void* flag, const int* flo, const int* fhi,
                                  const amrex_real* vfrac, const int* vlo, const int* vhi,
                                  const amrex_real* apx, const int* axlo, const int* axhi,
#if (AMREX_SPACEDIM >= 2)
                                  const amrex_real* apy, const int* aylo, const int* ayhi,
#if (AMREX_SPACEDIM == 3)
                                  const amrex_real* apz, const int* azlo, const int* azhi,
#endif
#endif
                                  const amrex_real* fcx, const int* cxlo, const int* cxhi,
#if (AMREX_SPACEDIM >= 2)
                                  const amrex_real* fcy, const int* cylo, const int* cyhi,
#if (AMREX_SPACEDIM == 3)
                                  const amrex_real* fcz, const int* czlo, const int* czhi,
#endif
#endif
                                  const amrex_real* dxinv,
                                  const amrex_real alpha, const amrex_real beta);

#ifdef __cplusplus
}
#endif

#endif
