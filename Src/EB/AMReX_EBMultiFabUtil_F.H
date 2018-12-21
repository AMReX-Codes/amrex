#ifndef AMREX_EBMULTIFABUTIL_F_H_
#define AMREX_EBMULTIFABUTIL_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void amrex_eb_set_covered (const int* lo, const int* hi,
                               amrex_real* data, const int* dlo, const int* dhi,
                               const void* flag, const int* flo, const int* fhi,
                               const amrex_real* v, const int* ncomp);

    void amrex_eb_set_covered_nodes (const int* lo, const int* hi,
                                     amrex_real* data, const int* dlo, const int* dhi,
                                     const void* flag, const int* flo, const int* fhi,
                                     const amrex_real* v, const int* ncomp);

    void amrex_eb_set_covered_faces (const int* lo, const int* hi,
                                     amrex_real* data, const int* dlo, const int* dhi,
                                     const amrex_real* a, const int* alo, const int* ahi);

    void amrex_eb_avgdown_sv (const int* lo, const int* hi,
                              const amrex_real* fine, const int* flo, const int* fhi,
                              amrex_real* crse, const int* clo, const int* chi,
                              const amrex_real* vol, const int* vlo, const int* vhi,
                              const amrex_real* vfol, const int* vflo, const int* vfhi,
                              const int* ratio, const int* ncomp);

    void amrex_eb_avgdown (const int* lo, const int* hi,
                           const amrex_real* fine, const int* flo, const int* fhi,
                           amrex_real* crse, const int* clo, const int* chi,
                           const amrex_real* vfol, const int* vflo, const int* vfhi,
                           const int* ratio, const int* ncomp);
   
    void amrex_eb_avgdown_faces (const int* lo, const int* hi, 
                                 const amrex_real* fine, const int* flo, const int* fhi, 
                                 amrex_real* crse, const int* clo, const int* chi,
                                 const amrex_real* ap, const int* aplo, const int* aphi, 
                                 const int* ratio, const int* idir, const int* ncomp);

    void amrex_eb_avgdown_boundaries (const int* lo, const int* hi, 
                                      const amrex_real* fine, const int* flo, const int* fhi, 
                                      amrex_real* crse, const int* clo, const int* chi,
                                      const amrex_real* ba, const int* blo, const int* bhi, 
                                      const int* ratio, const int* ncomp);

    void amrex_compute_eb_divergence (const int* lo, const int* hi,
                                   amrex_real* divu, const int* dlo, const int* dhi,
                                   const amrex_real* u, const int* ulo, const int* uhi,
#if (AMREX_SPACEDIM > 1)
                                   const amrex_real* v, const int* vlo, const int* vhi,
#endif
#if (AMREX_SPACEDIM == 3)
                                   const amrex_real* w, const int* wlo, const int* whi,
#endif
                                   const int* ccmask, const int* cmlo, const int* cmhi,
                                   const void* flag, const int* flo, const int* fhi,
                                   const amrex_real* vfrac, const int* klo, const int* khi,
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
                                   const amrex_real* dxinv);


    void amrex_eb_avg_fc_to_cc (const int* lo, const int* hi,
                                amrex_real* cc, const int* cclo, const int* cchi,
                                const amrex_real* fx, const int* fxlo, const int* fxhi,
#if (AMREX_SPACEDIM >= 2)
                                const amrex_real* fy, const int* fylo, const int* fyhi,
#if (AMREX_SPACEDIM == 3)
                                const amrex_real* fz, const int* fzlo, const int* fzhi,
#endif
#endif
                                const amrex_real* apx, const int* axlo, const int* axhi,
#if (AMREX_SPACEDIM >= 2)
                                const amrex_real* apy, const int* aylo, const int* ayhi,
#if (AMREX_SPACEDIM == 3)
                                const amrex_real* apz, const int* azlo, const int* azhi,
#endif
#endif
                                const void* flag, const int* flo, const int* fhi);
    
#ifdef __cplusplus
}
#endif

#endif
