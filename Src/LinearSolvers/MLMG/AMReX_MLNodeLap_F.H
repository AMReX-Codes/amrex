#ifndef AMREX_ML_NODE_LAP_F_H_
#define AMREX_ML_NODE_LAP_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

#if (AMREX_SPACEDIM == 2)
    void amrex_mlndlap_set_rz (const int* rz);
#endif

    int amrex_mlndlap_any_fine_sync_cells (const int* lo, const int* hi,
                                           const int* msk, const int* mlo, const int* mhi,
                                           const int* fine_flag);

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
