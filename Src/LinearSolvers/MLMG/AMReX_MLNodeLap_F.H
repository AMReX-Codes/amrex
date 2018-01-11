#ifndef AMREX_ML_NODE_LAP_F_H_
#define AMREX_ML_NODE_LAP_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

#if (AMREX_SPACEDIM == 2)
    void amrex_mlndlap_set_rz (const int* rz);
#endif

    void amrex_mlndlap_set_nodal_mask (const int* lo, const int* hi,
                                       int* nmsk, const int* nlo, const int* nhi,
                                       const int* cmsk, const int* clo, const int* chi);

    void amrex_mlndlap_set_dirichlet_mask (int* dmsk, const int* dlo, const int* dhi,
                                           const int* omsk, const int* olo, const int* ohi,
                                           const int* domlo, const int* domhi,
                                           const void* bclo, const void* bchi);

    void amrex_mlndlap_fixup_res_mask (const int* lo, const int* hi,
                                       int* rmsk, const int* rlo, const int* rhi,
                                       const int* cfmsk, const int* flo, const int* fhi);

    void amrex_mlndlap_set_dot_mask (const int* lo, const int* hi,
                                     amrex_real* dmsk, const int* dlo, const int* dhi,
                                     const int* omsk, const int* olo, const int* ohi,
                                     const int* domlo, const int* domhi,
                                     const void* bclo, const void* bchi);

    int amrex_mlndlap_any_fine_sync_cells (const int* lo, const int* hi,
                                           const int* msk, const int* mlo, const int* mhi,
                                           const int* fine_flag);

    void amrex_mlndlap_avgdown_coeff (const int* lo, const int* hi,
                                      amrex_real* crse, const int* clo, const int* chi,
                                      const amrex_real* fine, const int* flo, const int* fhi,
                                      const int* idim);

    void amrex_mlndlap_fillbc_cc (amrex_real* sigma, const int* slo, const int* shi,
                                  const int* dlo, const int* dhi,
                                  const void* bclo, const void* bchi);

    void amrex_mlndlap_fillbc_cc_i (int* sigma, const int* slo, const int* shi,
                                    const int* dlo, const int* dhi,
                                    const void* bclo, const void* bchi);

    void amrex_mlndlap_applybc (amrex_real* phi, const int* hlo, const int* hhi,
                                const int* dlo, const int* dhi,
                                const void* bclo, const void* bchi);

    void amrex_mlndlap_impose_neumann_bc (const int* lo, const int* hi,
                                          amrex_real* rhs, const int* rlo, const int* rhi,
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
                              const int* dmsk, const int* dmlo, const int* dmhi,
                              const amrex_real* dxinv, const int* domlo, const int* domhi,
                              const void* bclo, const void* bchi);


    void amrex_mlndlap_adotx_aa (const int* lo, const int* hi,
                              amrex_real* y, const int* ylo, const int* yhi,
                              const amrex_real* x, const int* xlo, const int* xhi,
                              const amrex_real* sig, const int* slo, const int* shi,
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

    void amrex_mlndlap_divu (const int* lo, const int* hi,
                             amrex_real* rhs, const int* rlo, const int* rhi,
                             const amrex_real* vel, const int* vlo, const int* vhi,
                             const int* dmsk, const int* dmlo, const int* dmhi,
                             const amrex_real* dxinv);

    void amrex_mlndlap_rhcc (const int* lo, const int* hi,
                             amrex_real* rhs, const int* rlo, const int* rhi,
                             const amrex_real* rhcc, const int* clo, const int* chi,
                             const int* dmsk, const int* mlo, const int* mhi);

    void amrex_mlndlap_mknewu (const int* lo, const int* hi,
                               amrex_real* u, const int* ulo, const int* uhi,
                               const amrex_real* p, const int* plo, const int* phi,
                               const amrex_real* sigma, const int* slo, const int* shi,
                               const amrex_real* dxinv);

    void amrex_mlndlap_divu_fine_contrib (const int* clo, const int* chi,
                                          const int* cgrid_lo, const int* cgrid_hi,
                                          amrex_real* rhs, const int* rlo, const int* rhi,
                                          const amrex_real* vel, const int* vlo, const int* vhi,
                                          amrex_real* frh, const int* flo, const int* fhi,
                                          const int* dmsk, const int* mlo, const int* mhi,
                                          const amrex_real* dxinv);

    void amrex_mlndlap_divu_cf_contrib (const int* lo, const int* hi,
                                        amrex_real* rhs, const int* rlo, const int* rhi,
                                        const amrex_real* vel, const int* vlo, const int* vhi,
                                        const int* dmsk, const int* dmlo, const int* dmhi,
                                        const int* ndmsk, const int* ndmlo, const int* ndmhi,
                                        const int* ccmsk, const int* ccmlo, const int* ccmhi,
                                        const amrex_real* cc, const int* clo, const int* chi,
                                        const amrex_real* dxinv, const int* ndlo, const int* ndhi,
                                        const void* bclo, const void* bchi);

    void amrex_mlndlap_crse_resid (const int* lo, const int* hi,
                                   amrex_real* resid, const int* rslo, const int* rshi,
                                   const amrex_real* rhs, const int* rhlo, const int* rhhi,
                                   const int* msk, const int* mlo, const int* mhi,
                                   const int* ndlo, const int* ndhi,
                                   const void* bclo, const void* bchi);

    void amrex_mlndlap_res_fine_contrib (const int* clo, const int* chi,
                                         const int* cgrid_lo, const int* cgrid_hi,
                                         amrex_real* fc, const int* flo, const int* fhi,
                                         const amrex_real* phi, const int* phlo, const int* phhi,
                                         const amrex_real* sig, const int* slo, const int* shi,
                                         amrex_real* Ax, const int* alo, const int* ahi,
                                         const int* dmsk, const int* dmlo, const int* dmhi,
                                         const amrex_real* dxinv);

    void amrex_mlndlap_res_cf_contrib (const int* lo, const int* hi,
                                       amrex_real* res, const int* rlo, const int* rhi,
                                       const amrex_real* phi, const int* phlo, const int* phhi,
                                       const amrex_real* rhs, const int* rhlo, const int* rhhi,
                                       const amrex_real* sig, const int* slo, const int* shi,
                                       const int* dmsk, const int* mlo, const int* mhi,
                                       const int* ndmsk, const int* ndmlo, const int* ndmhi,
                                       const int* ccmsk, const int* ccmlo, const int* ccmhi,
                                       const amrex_real* cc, const int* clo, const int* chi,
                                       const amrex_real* dxinv, const int* ndlo, const int* ndhi,
                                       const void* bclo, const void* bchi);

    void amrex_mlndlap_zero_fine (const int* lo, const int* hi,
                                  amrex_real* phi, const int* philo, const int* phihi,
                                  const int* msk, const int* mlo, const int* mhi,
                                  const int* fine_flag);


    void amrex_mlndlap_rhcc_fine_contrib (const int* clo, const int* chi,
                                          const int* cgrid_lo, const int* cgrid_hi,
                                          amrex_real* rhs, const int* rlo, const int* rhi,
                                          const amrex_real* cc, const int* cclo, const int* cchi,
                                          const int* dmsk, const int* mlo, const int* mhi);

    void amrex_mlndlap_rhcc_crse_contrib (const int* lo, const int* hi,
                                          amrex_real* crhs, const int* rlo, const int* rhi,
                                          const amrex_real* cc, const int* cclo, const int* cchi,
                                          const int* dmsk, const int* dmlo, const int* dmhi,
                                          const int* ndmsk, const int* ndmlo, const int* ndmhi,
                                          const int* ccmsk, const int* ccmlo, const int* ccmhi);

#ifdef __cplusplus
}
#endif

#endif
