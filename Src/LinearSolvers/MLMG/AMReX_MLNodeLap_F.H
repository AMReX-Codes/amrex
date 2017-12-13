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
                                     const int* dlo, const int* dhi,
                                     const void* bclo, const void* bchi);

    void amrex_mlndlap_divu (const int* lo, const int* hi,
                             amrex_real* rhs, const int* rlo, const int* rhi,
                             const amrex_real* vel, const int* vlo, const int* vhi,
                             const int* dmsk, const int* dmlo, const int* dmhi,
                             const amrex_real* dxinv, const int* ndlo, const int* ndhi,
                             const void* bclo, const void* bchi);

    void amrex_mlndlap_mknewu (const int* lo, const int* hi,
                               amrex_real* u, const int* ulo, const int* uhi,
                               const amrex_real* p, const int* plo, const int* phi,
                               const amrex_real* sigma, const int* slo, const int* shi,
                               const amrex_real* dxinv);

    void amrex_mlndlap_applybc (amrex_real* phi, const int* hlo, const int* hhi,
                                const int* dlo, const int* dhi,
                                const void* bclo, const void* bchi);

    void amrex_mlndlap_adotx_ha (const int* lo, const int* hi,
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
                              const int* dmsk, const int* dmlo, const int* dmhi,
                              const amrex_real* dxinv, const int* domlo, const int* domhi,
                              const void* bclo, const void* bchi);


    void amrex_mlndlap_adotx_aa (const int* lo, const int* hi,
                              amrex_real* y, const int* ylo, const int* yhi,
                              const amrex_real* x, const int* xlo, const int* xhi,
                              const amrex_real* sig, const int* slo, const int* shi,
                              amrex_real* dg, const int* dlo, const int* dhi,
                              const int* dmsk, const int* dmlo, const int* dmhi,
                              const amrex_real* dxinv, const int* domlo, const int* domhi,
                              const void* bclo, const void* bchi);


    void amrex_mlndlap_gauss_seidel_ha (const int* lo, const int* hi,
                                     amrex_real* sol, const int* slo, const int* shi,
                                     const amrex_real* rhs, const int* rlo, const int* rhi,
                                     const amrex_real* sx, const int* sxlo, const int* sxhi,
#if (AMREX_SPACEDIM >= 2)
                                     const amrex_real* sy, const int* sylo, const int* syhi,
#if (AMREX_SPACEDIM == 3)
                                     const amrex_real* sz, const int* szlo, const int* szhi,
#endif
#endif
                                     const int* dmsk, const int* dmlo, const int* dmhi,
                                     const amrex_real* dxinv, const int* domlo, const int* domhi,
                                     const void* bclo, const void* bchi);


    void amrex_mlndlap_gauss_seidel_aa (const int* lo, const int* hi,
                                     amrex_real* sol, const int* slo, const int* shi,
                                     const amrex_real* rhs, const int* rlo, const int* rhi,
                                     const amrex_real* sig, const int* sglo, const int* sghi,
                                     const int* dmsk, const int* dmlo, const int* dmhi,
                                     const amrex_real* dxinv, const int* domlo, const int* domhi,
                                     const void* bclo, const void* bchi);


    void amrex_mlndlap_jacobi_ha (const int* lo, const int* hi,
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
                               const int* dmsk, const int* dmlo, const int* dmhi,
                               const amrex_real* dxinv, const int* domlo, const int* domhi,
                               const void* bclo, const void* bchi);
                               

    void amrex_mlndlap_jacobi_aa (const int* lo, const int* hi,
                               amrex_real* sol, const int* slo, const int* shi,
                               const amrex_real* Ax, const int* alo, const int* ahi,
                               const amrex_real* rhs, const int* rlo, const int* rhi,
                               const amrex_real* sig, const int* sglo, const int* sghi,
                               const int* dmsk, const int* dmlo, const int* dmhi,
                               const amrex_real* dxinv, const int* domlo, const int* domhi,
                               const void* bclo, const void* bchi);
                               

    void amrex_mlndlap_restriction (const int* lo, const int* hi,
                                    amrex_real* crse, const int* clo, const int* chi,
                                    const amrex_real* fine, const int* flo, const int* fhi,
                                    const int* dmsk, const int* dmlo, const int* dmhi,
                                    const int* dlo, const int* dhi,
                                    const void* bclo, const void* bchi);

    void amrex_mlndlap_interpolation_ha (const int* clo, const int* chi,
                                      amrex_real* fine, const int* fflo, const int* ffhi,
                                      const amrex_real* crse, const int* cflo, const int* cfhi,
                                      const amrex_real* sx, const int* sxlo, const int* sxhi,
#if (AMREX_SPACEDIM >= 2)
                                      const amrex_real* sy, const int* sylo, const int* syhi,
#if (AMREX_SPACEDIM == 3)
                                      const amrex_real* sz, const int* szlo, const int* szhi,
#endif
#endif
                                      const int* dmsk, const int* dmlo, const int* dmhi,
                                      const int* domlo, const int* domhi,
                                      const void* bclo, const void* bchi);


    void amrex_mlndlap_interpolation_aa (const int* clo, const int* chi,
                                      amrex_real* fine, const int* fflo, const int* ffhi,
                                      const amrex_real* crse, const int* cflo, const int* cfhi,
                                      const amrex_real* sig, const int* siglo, const int* sighi,
                                      const int* dmsk, const int* dmlo, const int* dmhi,
                                      const int* domlo, const int* domhi,
                                      const void* bclo, const void* bchi);


    void amrex_mlndlap_zero_fine (const int* lo, const int* hi,
                                  amrex_real* phi, const int* philo, const int* phihi,
                                  const int* msk, const int* mlo, const int* mhi);

    void amrex_mlndlap_crse_resid (const int* lo, const int* hi,
                                   amrex_real* resid, const int* rslo, const int* rshi,
                                   const amrex_real* rhs, const int* rhlo, const int* rhhi,
                                   const int* msk, const int* mlo, const int* mhi);

    int amrex_mlndlap_any_zero (const int* lo, const int* hi,
                                const int* msk, const int* mlo, const int* mhi);

    void amrex_mlndlap_set_dirichlet_mask (int* dmsk, const int* dlo, const int* dhi,
                                           const int* omsk, const int* olo, const int* ohi,
                                           const int* domlo, const int* domhi,
                                           const void* bclo, const void* bchi);

#ifdef __cplusplus
}
#endif

#endif
