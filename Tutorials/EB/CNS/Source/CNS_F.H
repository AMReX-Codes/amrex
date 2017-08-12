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

#ifdef __cplusplus
}
#endif

#endif
