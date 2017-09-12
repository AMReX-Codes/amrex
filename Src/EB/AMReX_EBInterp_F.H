#ifndef AMREX_EBINTERP_F_H_
#define AMREX_EBINTERP_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_ebinterp_pc_sv (const int* tflo, const int* tfhi,
                               const int* tclo, const int* tchi,
                               const amrex_real* crse, const int* clo, const int* chi,
                               amrex_real* fine, const int* flo, const int* fhi,
                               const int* ncomp, const int* ratio,
                               const int* cdlo, const int* cdhi,
                               const void* flag, const int* fglo, const int* fghi);

#ifdef __cplusplus
}
#endif

#endif
