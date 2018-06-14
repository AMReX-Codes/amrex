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

    void amrex_mlebabeclap_gsrb  (const int* lo, const int* hi,
                                  amrex_real* sol, const int* slo, const int* shi,
                                  const amrex_real* rhs, const int* rlo, const int* rhi,
                                  const amrex_real* a, const int* alo, const int* ahi,
                                  const amrex_real* bx, const int* bxlo, const int* bxhi,
#if (AMREX_SPACEDIM >= 2)
                                  const amrex_real* by, const int* bylo, const int* byhi,
#if (AMREX_SPACEDIM == 3)
                                  const amrex_real* bz, const int* bzlo, const int* bzhi,
#endif
#endif
                                  const int* m0, const int* m0lo, const int* m0hi,
#if (AMREX_SPACEDIM >= 2)
                                  const int* m2, const int* m2lo, const int* m2hi,
#if (AMREX_SPACEDIM == 3)
                                  const int* m4, const int* m4lo, const int* m4hi,
#endif
#endif
                                  const int* m1, const int* m1lo, const int* m1hi,
#if (AMREX_SPACEDIM >= 2)
                                  const int* m3, const int* m3lo, const int* m3hi,
#if (AMREX_SPACEDIM == 3)
                                  const int* m5, const int* m5lo, const int* m5hi,
#endif
#endif
                                  const amrex_real* f0, const int* f0lo, const int* f0hi,
#if (AMREX_SPACEDIM >= 2)
                                  const amrex_real* f2, const int* f2lo, const int* f2hi,
#if (AMREX_SPACEDIM == 3)
                                  const amrex_real* f4, const int* f4lo, const int* f4hi,
#endif
#endif
                                  const amrex_real* f1, const int* f1lo, const int* f1hi,
#if (AMREX_SPACEDIM >= 2)
                                  const amrex_real* f3, const int* f3lo, const int* f3hi,
#if (AMREX_SPACEDIM == 3)
                                  const amrex_real* f5, const int* f5lo, const int* f5hi,
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
                                  const amrex_real alpha, const amrex_real beta, const int redblack);


    void amrex_mlebabeclap_normalize (const int* lo, const int* hi,
                                      amrex_real* y, const int* ylo, const int* yhi,
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


    void amrex_eb_mg_interp (const int* lo, const int* hi,
                             amrex_real* fine, const int* flo, const int* fhi,
                             const amrex_real* crse, const int* clo, const int* chi,
                             const void* flag, const int* glo, const int* ghi,
                             const int* ncomp);

#ifdef __cplusplus
}
#endif

#endif
