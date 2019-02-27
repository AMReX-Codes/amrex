#ifndef AMREX_MLMG_F_H_
#define AMREX_MLMG_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_mlmg_lin_nd_interp (const int* clo, const int* chi,
                                   const int* flo, const int* fhi,
                                   amrex_real* fine, const int* fdlo, const int* fdhi,
                                   const amrex_real* crse, const int* cdlo, const int* cdhi,
				   const int* nc);

#ifdef AMREX_USE_EB
    void amrex_mlmg_eb_cc_interp (const int* lo, const int* hi,
                                  amrex_real* fine, const int* fdlo, const int* fdhi,
                                  const amrex_real* crse, const int* cdlo, const int* cdhi,
                                  const void* flag, const int* glo, const int* ghi,
                                  const int* ratio,const int* nc);
#endif

#ifdef __cplusplus
}
#endif

#endif
