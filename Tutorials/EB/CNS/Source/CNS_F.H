#ifndef CNS_F_H_
#define CNS_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void cns_init_fort (const int* phys_bc_lo, const int* phys_bc_hi,
                        int, int, int, int, int, int,
                        const amrex_real* problo, const amrex_real* probhi);

    void cns_hypfill (BL_FORT_FAB_ARG_3D(state),
                      const int* dlo, const int* dhi,
                      const amrex::Real* dx, const amrex::Real* glo, 
                      const amrex::Real* time, const int* bc);
    
    void cns_denfill (BL_FORT_FAB_ARG_3D(state),
                      const int* dlo, const int* dhi,
                      const amrex::Real* dx, const amrex::Real* glo, 
                      const amrex::Real* time, const int* bc);

    void cns_derpres (BL_FORT_FAB_ARG_3D(der),const int* nvar,
                      const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
                      const int* lo, const int* hi,
                      const int* domain_lo, const int* domain_hi,
                      const Real* delta, const Real* xlo,
                      const Real* time, const Real* dt, const int* bcrec, 
                      const int* level, const int* grid_no);

    void cns_dervel (BL_FORT_FAB_ARG_3D(der),const int* nvar,
                     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
                     const int* lo, const int* hi,
                     const int* domain_lo, const int* domain_hi,
                     const Real* delta, const Real* xlo,
                     const Real* time, const Real* dt, const int* bcrec, 
                     const int* level, const int* grid_no);
    
#ifdef __cplusplus
}
#endif

#endif
