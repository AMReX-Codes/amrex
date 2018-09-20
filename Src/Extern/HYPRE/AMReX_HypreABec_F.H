#ifndef AMREX_HYPREABEC_F_H_
#define AMREX_HYPREABEC_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_hpacoef (const int* lo, const int* hi, amrex_real* mat,
                        const amrex_real* acoefs, const int* alo, const int* ahi,
                        const amrex_real* sa);

    void amrex_hpbcoef (const int* lo, const int* hi, amrex_real* mat,
                        const amrex_real* bcoefs, const int* blo, const int* bhi,
                        const amrex_real* sb, const amrex_real* dx, const int* idim);

    void amrex_hpmat (const int* lo, const int* hi, amrex_real* mat,
                      const amrex_real* bcoefs, const int* blo, const int* bhi,
                      const int* mask, const int* mlo, const int* mhi,
                      const amrex_real* sb, const amrex_real* dx, const int* cdir,
                      const int* bct, const amrex_real* bcl, const int* bho);

    void amrex_hpdiag (const int* lo, const int* hi, amrex_real* mat,
                       amrex_real* diag, const int* dlo, const int* dhi);

    void amrex_hpijmatrix (const int* lo, const int* hi,
                           const void* nrows, void* ncols,
                           void* rows, void* cols, amrex_real* mat,
                           const void* cell_id, const int* clo, const int* chi,
                           const void* cell_id_begin,
                           amrex_real* diaginv, const int* dlo, const int* dhi,
                           const amrex_real* ac, const int* alo, const int* ahi,
                           const amrex_real* bx, const int* bxlo, const int* bxhi,
#if (AMREX_SPACEDIM >= 2)
                           const amrex_real* by, const int* bylo, const int* byhi,
#if (AMREX_SPACEDIM == 3)
                           const amrex_real* bz, const int* bzlo, const int* bzhi,
#endif
#endif
                           const amrex_real* sa, const amrex_real* sb, const amrex_real* dx,
                           const int* bct, const amrex_real* bcl, const int* bho);
                           

#ifdef AMREX_USE_EB
    void amrex_hpeb_fill_cellid (const int* lo, const int* hi,
                                 void* nrows, const void* cell_id, const int* clo, const int* chi,
                                 const void* flag, const int* flo, const int* fhi);

    void amrex_hpeb_copy_from_vec (const int* lo, const int* hi,
                                   amrex_real* a, const int* alo, const int* ahi,
                                   const amrex_real* v, const void* nv,
                                   const void* flag, const int* flo, const int* fhi);

    void amrex_hpeb_copy_to_vec (const int* lo, const int* hi,
                                 const amrex_real* a, const int* alo, const int* ahi,
                                 amrex_real* v, const void* nv,
                                 const void* flag, const int* flo, const int* fhi);

    void amrex_hpeb_ijmatrix (const int* lo, const int* hi,
                              const void* nrows, void* ncols,
                              void* rows, void* cols, amrex_real* mat,
                              const void* cell_id, const int* clo, const int* chi,
                              const void* cell_id_begin,
                              amrex_real* diaginv, const int* dlo, const int* dhi,
                              const amrex_real* ac, const int* alo, const int* ahi,
                              const amrex_real* bx, const int* bxlo, const int* bxhi,
#if (AMREX_SPACEDIM >= 2)
                              const amrex_real* by, const int* bylo, const int* byhi,
#if (AMREX_SPACEDIM == 3)
                              const amrex_real* bz, const int* bzlo, const int* bzhi,
#endif
#endif
                              const void* flag, const int* flo, const int* fhi,
                              const amrex_real* vfrc, const int* vlo, const int* vhi,
                              const amrex_real* apx, const int* axlo, const int* axhi,
#if (AMREX_SPACEDIM >= 2)
                              const amrex_real* apy, const int* aylo, const int* ayhi,
#if (AMREX_SPACEDIM == 3)
                              const amrex_real* apz, const int* azlo, const int* azhi,
#endif
#endif
                              const amrex_real* fcx, const int* fxlo, const int* fxhi,
#if (AMREX_SPACEDIM >= 2)
                              const amrex_real* fcy, const int* fylo, const int* fyhi,
#if (AMREX_SPACEDIM == 3)
                              const amrex_real* fcz, const int* fzlo, const int* fzhi,
#endif
#endif
                              const amrex_real* ba,   const int* balo, const int* bahi,
                              const amrex_real* bcen, const int* bclo, const int* bchi,
                              const amrex_real* beb,  const int* elo,  const int* ehi,
                              const int* is_eb_dirichlet,
                              const amrex_real* sa, const amrex_real* sb, const amrex_real* dx,
                              const int* bct, const amrex_real* bcl, const int* bho);

#endif
    
#ifdef __cplusplus
};
#endif

#endif /* _HYPREABEC_F_H_ */
