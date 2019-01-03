#ifndef AMREX_FLUXREG_F_H_
#define AMREX_FLUXREG_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void amrex_froverwrite_cfb (const int* lo, const int* hi,
                                amrex_real* dst, const int* dlo, const int* dhi,
                                const amrex_real* src, const int* slo, const int* shi,
                                const int* mask, const int* mlo, const int* mhi,
                                const int* ncomp, const int* idim, const amrex_real* scale);

#ifdef __cplusplus
}
#endif

#endif
