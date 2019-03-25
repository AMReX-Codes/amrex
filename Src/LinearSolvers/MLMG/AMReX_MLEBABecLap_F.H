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
                                  const int* ccmask, const int* cmlo, const int* cmhi,
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
                                  const amrex_real* ba, const int* balo, const int* bahi,
                                  const amrex_real* bc, const int* bclo, const int* bchi,
                                  const amrex_real* beb, const int* elo, const int* ehi,
                                  const int is_dirichlet, const int is_ho_dirichlet,
                                  const amrex_real* phieb, const int* plo, const int* phi,
                                  const int is_inhomog, const amrex_real* dxinv,
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
                                  const int* ccmask, const int* cmlo, const int* cmhi,
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
                                  const amrex_real* ba, const int* balo, const int* bahi,
                                  const amrex_real* bc, const int* bclo, const int* bchi,
                                  const amrex_real* beb, const int* elo, const int* ehi,
                                  const int is_dirichlet, const int is_ho_dirichlet,
                                  const amrex_real* dxinv,
                                  const amrex_real alpha, const amrex_real beta, const int redblack);

    void amrex_mlebabeclap_flux( const int* lo, const int* hi, 
                                 amrex_real* fx, const int* fxlo, const int* fxhi, 
#if (AMREX_SPACEDIM >= 2) 
                                 amrex_real* fy, const int* fylo, const int* fyhi, 
#if (AMREX_SPACEDIM == 3) 
                                 amrex_real* fz, const int* fzlo, const int* fzhi, 
#endif
#endif
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

                                 const amrex_real* sol, const int* slo, const int* shi,
                                 const amrex_real* bx, const int* bxlo, const int* bxhi,
#if (AMREX_SPACEDIM >= 2)
                                 const amrex_real* by, const int* bylo, const int* byhi,
#if (AMREX_SPACEDIM == 3)
                                 const amrex_real* bz, const int* bzlo, const int* bzhi,
#endif
#endif
                                 const int* ccmask, const int* cmlo, const int* cmhi,
                                 const void* flag, const int* flo, const int* fhi,
                                 const amrex_real* dxinv, const amrex_real beta,
                                 const int face_only);

    void amrex_mlebabeclap_grad( const int* xlo, const int* xhi, 
#if (AMREX_SPACEDIM >= 2) 
                                 const int* ylo, const int* yhi, 
#if (AMREX_SPACEDIM == 3) 
                                 const int* zlo, const int* zhi, 
#endif
#endif // */
                                 const amrex_real* dat, const int* dlo, const int* dhi,
                                 amrex_real* gx, const int* gxlo, const int* gxhi,
#if (AMREX_SPACEDIM >= 2)
                                 amrex_real* gy, const int* gylo, const int* gyhi,
#if (AMREX_SPACEDIM == 3)
                                 amrex_real* gz, const int* gzlo, const int* gzhi,
#endif
#endif
                                 const amrex_real* ax, const int* axlo, const int* axhi, 
#if (AMREX_SPACEDIM >= 2)
                                 const amrex_real* ay, const int* aylo, const int* ayhi, 
#if (AMREX_SPACEDIM == 3)
                                 const amrex_real* az, const int* azlo, const int* azhi,
#endif
#endif
                                 const amrex_real* fcx, const int* fcxlo, const int* fcxhi, 
#if (AMREX_SPACEDIM >= 2)
                                 const amrex_real* fcy, const int* fcylo, const int* fcyhi, 
#if (AMREX_SPACEDIM == 3)
                                 const amrex_real* fcz, const int* fczlo, const int* fczhi,
#endif
#endif
                                 const int* ccmask, const int* cmlo, const int* cmhi,
                                 const void* flag, const int* flo, const int* fhi,
                                 const amrex_real* dxinv);

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
                                      const int* ccmask, const int* cmlo, const int* cmhi,
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
                                      const amrex_real* ba, const int* balo, const int* bahi,
                                      const amrex_real* bc, const int* bclo, const int* bchi,
                                      const amrex_real* beb, const int* elo, const int* ehi,
                                      const int is_dirichlet, const int is_ho_dirichlet,
                                      const amrex_real* dxinv,
                                      const amrex_real alpha, const amrex_real beta);


    void amrex_eb_mg_interp (const int* lo, const int* hi,
                             amrex_real* fine, const int* flo, const int* fhi,
                             const amrex_real* crse, const int* clo, const int* chi,
                             const void* flag, const int* glo, const int* ghi,
                             const int* ncomp);

    void amrex_mlebabeclap_apply_bc (const int* lo, const int* hi,
                                     amrex_real* phi, const int* philo, const int* phihi,
                                     const void* flag, const int* flo, const int* fhi,
                                     const amrex_real* apx, const int* axlo, const int* axhi,
#if (AMREX_SPACEDIM >= 2)
                                     const amrex_real* apy, const int* aylo, const int* ayhi,
#if (AMREX_SPACEDIM == 3)
                                     const amrex_real* apz, const int* azlo, const int* azhi,
#endif
#endif
                                     const int* mask, const int* mlo, const int* mhi,
                                     const int cdir, const int bct, const amrex_real bcl,
                                     const amrex_real* bcval, const int* blo, const int* bhi,
                                     const int maxorder, const amrex_real* dxinv,
                                     const int inhomog, const int nc);

    void amrex_eb_copy_dirichlet (const int* lo, const int* hi,
                                  amrex_real* phi, const int* phlo, const int* phhi,
                                  const amrex_real* phiin, const int* phinlo, const int* phinhi,
                                  amrex_real* beta, const int* blo, const int* bhi,
                                  const amrex_real* betain, const int* binlo, const int* binhi,
                                  const void* flag, const int* flo, const int* fhi);

    void amrex_eb_homog_dirichlet (const int* lo, const int* hi,
                                   amrex_real* phi, const int* phlo, const int* phhi,
                                   amrex_real* beta, const int* blo, const int* bhi,
                                   const amrex_real* betain, const int* binlo, const int* binhi,
                                   const void* flag, const int* flo, const int* fhi);

#ifdef __cplusplus
}
#endif

#endif
