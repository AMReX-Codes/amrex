#ifndef AMREX_ML_NODE_LAP_F_H_
#define AMREX_ML_NODE_LAP_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_mlndlap_avgdown_coeff (const int* lo, const int* hi,
                                      amrex_real* crse, const int* clo, const int* chi,
                                      const amrex_real* fine, const int* flo, const int* fhi,
                                      const int* idim);

    void amrex_mlndlap_fillbc_coeff (amrex_real* sigma, const int* slo, const int* shi,
                                     const int* dlo, const int* dhi);

    void amrex_mlndlap_divu (const int* lo, const int* hi,
                             amrex_real* rhs, const int* rlo, const int* rhi,
                             const amrex_real* vel, const int* vlo, const int* vhi,
                             const amrex_real* dxinv, const int* ndlo, const int* ndhi,
                             const void* bclo, const void* bchi);

    void amrex_mlndlap_applybc (amrex_real* phi, const int* hlo, const int* hhi,
                                const int* dlo, const int* dhi,
                                const void* bclo, const void* bchi);

    void amrex_mlndlap_adotx (const int* lo, const int* hi,
                              amrex_real* y, const int* ylo, const int* yhi,
                              const amrex_real* x, const int* xlo, const int* xhi,
                              const amrex_real* sx, const int* sxlo, const int* sxhi,
#if (AMREX_SPACEDIM >= 2)
                              const amrex_real* sy, const int* sylo, const int* syhi,
#if (AMREX_SPACEDIM == 3)
                              const amrex_real* sz, const int* szlo, const int* szhi,
#endif
#endif
                              amrex_real* dg, const int* dlo, const int* dhi,
                              const amrex_real* dxinv);


    void amrex_mlndlap_jacobi (const int* lo, const int* hi,
                               amrex_real* sol, const int* slo, const int* shi,
                               const amrex_real* Ax, const int* alo, const int* ahi,
                               const amrex_real* rhs, const int* rlo, const int* rhi,
                               const amrex_real* sx, const int* sxlo, const int* sxhi,
#if (AMREX_SPACEDIM >= 2)
                               const amrex_real* sy, const int* sylo, const int* syhi,
#if (AMREX_SPACEDIM == 3)
                               const amrex_real* sz, const int* szlo, const int* szhi,
#endif
#endif
                               const amrex_real* dxinv);
                               

    void amrex_mlndlap_restriction (const int* lo, const int* hi,
                                    amrex_real* crse, const int* clo, const int* chi,
                                    const amrex_real* fine, const int* flo, const int* fhi,
                                    const int* dlo, const int* dhi,
                                    const void* bclo, const void* bchi);

    void amrex_mlndlap_interpolation (const int* clo, const int* chi,
                                      amrex_real* fine, const int* fflo, const int* ffhi,
                                      const amrex_real* crse, const int* cflo, const int* cfhi,
                                      const amrex_real* sx, const int* sxlo, const int* sxhi
#if (AMREX_SPACEDIM >= 2)
                                      , const amrex_real* sy, const int* sylo, const int* syhi
#if (AMREX_SPACEDIM == 3)
                                      , const amrex_real* sz, const int* szlo, const int* szhi
#endif
#endif
        );

#ifdef __cplusplus
}
#endif

#endif
