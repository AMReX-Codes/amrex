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

    void amrex_eb2_build_faces (const int* lo, const int* hi,
                                void* cell, const int* clo, const int* chi,
                                void* fx, const int* fxlo, const int* fxhi,
                                void* fy, const int* fylo, const int* fyhi,
                                void* fz, const int* fzlo, const int* fzhi,
                                void* ex, const int* exlo, const int* exhi,
                                void* ey, const int* eylo, const int* eyhi,
                                void* ez, const int* ezlo, const int* ezhi,
                                const amrex_real* s, const int* slo, const int* shi,
                                const amrex_real* interx, const int* ixlo, const int* ixhi,
                                const amrex_real* intery, const int* iylo, const int* iyhi,
                                const amrex_real* interz, const int* izlo, const int* izhi,
                                amrex_real* apx, const int* axlo, const int* axhi,
                                amrex_real* apy, const int* aylo, const int* ayhi,
                                amrex_real* apz, const int* azlo, const int* azhi,
                                amrex_real* fcx, const int* cxlo, const int* cxhi,
                                amrex_real* fcy, const int* cylo, const int* cyhi,
                                amrex_real* fcz, const int* czlo, const int* czhi,
                                const amrex_real* dx, const amrex_real* dxinv, const amrex_real* problo);

    void amrex_eb2_build_cells (const int* lo, const int* hi,
                                void* cell, const int* clo, const int* chi,
                                void* fx, const int* fxlo, const int* fxhi,
                                void* fy, const int* fylo, const int* fyhi,
                                void* fz, const int* fzlo, const int* fzhi,
                                amrex_real* apx, const int* axlo, const int* axhi,
                                amrex_real* apy, const int* aylo, const int* ayhi,
                                amrex_real* apz, const int* azlo, const int* azhi,
                                const amrex_real* fcx, const int* cxlo, const int* cxhi,
                                const amrex_real* fcy, const int* cylo, const int* cyhi,
                                const amrex_real* fcz, const int* czlo, const int* czhi,
                                amrex_real* vfrac, const int* vlo, const int* vhi,
                                amrex_real* vcent, const int* tlo, const int* thi,
                                amrex_real* barea, const int* alo, const int* ahi,
                                amrex_real* bcent, const int* blo, const int* bhi,
                                amrex_real* bnorm, const int* mlo, const int* mhi);

#endif

#ifdef __cplusplus
}
#endif

#endif
