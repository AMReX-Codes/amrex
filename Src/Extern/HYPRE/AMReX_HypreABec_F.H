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

    void amrex_hpijmatrix (const int* lo, const int* hi,
                           const void* nrows, void* ncols,
                           void* rows, void* cols, amrex_real* mat,
                           const void* cell_id, const int* clo, const int* chi,
                           const void* cell_id_begin,
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
                           
    
#ifdef __cplusplus
};
#endif

#endif /* _HYPREABEC_F_H_ */
