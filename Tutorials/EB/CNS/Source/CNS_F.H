#ifndef CNS_F_H_
#define CNS_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void cns_init_fort (amrex_real gamma, amrex_real mu, amrex_real Pr, amrex_real Sc,
                        amrex_real C_S, amrex_real T_S,
                        const int* phys_bc_lo, const int* phys_bc_hi,
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
    
#ifdef __cplusplus
}
#endif

#endif
