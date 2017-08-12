#ifndef CNS_F_H_
#define CNS_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void cns_initdata (const int* level, const amrex_real* time,
                       const int* lo, const int* hi,
                       amrex_real* u, const int* ulo, const int* uhi,
                       const amrex_real* dx, const amrex_real* problo);

    void cns_init_fort (const int* phys_bc_lo, const int* phys_bc_hi,
                        int, int, int, int, int, int,
                        const amrex_real* problo, const amrex_real* probhi);

    void cns_hypfill (BL_FORT_FAB_ARG_3D(state),
                      const int* dlo, const int* dhi,
                      const amrex_real* dx, const amrex_real* glo, 
                      const amrex_real* time, const int* bc);
    
    void cns_denfill (BL_FORT_FAB_ARG_3D(state),
                      const int* dlo, const int* dhi,
                      const amrex_real* dx, const amrex_real* glo, 
                      const amrex_real* time, const int* bc);

    void cns_derpres (BL_FORT_FAB_ARG_3D(der),const int* nvar,
                      const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
                      const int* lo, const int* hi,
                      const int* domain_lo, const int* domain_hi,
                      const amrex_real* delta, const amrex_real* xlo,
                      const amrex_real* time, const amrex_real* dt, const int* bcrec, 
                      const int* level, const int* grid_no);

    void cns_dervel (BL_FORT_FAB_ARG_3D(der),const int* nvar,
                     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
                     const int* lo, const int* hi,
                     const int* domain_lo, const int* domain_hi,
                     const amrex_real* delta, const amrex_real* xlo,
                     const amrex_real* time, const amrex_real* dt, const int* bcrec, 
                     const int* level, const int* grid_no);

    void cns_compute_temperature (const int* lo, const int* hi,
                                  amrex_real* u, const int* ulo, const int* uhi);

    void cns_estdt (const int* lo, const int* hi,
                    const amrex_real* u, const int* ulo, const int* uhi,
                    const amrex_real* dx, amrex_real* dt);
                    

#ifdef __cplusplus
}
#endif

#endif
