#ifndef AMREX_LO_F_H_
#define AMREX_LO_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void amrex_lo_applybc(
        const int *flagden, const int *flagbc, const int *maxorder,
        amrex_real *phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
        const int *cdr,
        const int *bct,
        const amrex_real *bcl,
        const amrex_real *bcval, ARLIM_P(bcval_lo), ARLIM_P(bcval_hi),
        const int *mask,   ARLIM_P(mask_lo),  ARLIM_P(mask_hi), 
        amrex_real *den,         ARLIM_P(den_lo),   ARLIM_P(den_hi),
        const int *lo, const int *hi, const int *nc,
        const amrex_real *h
        );

    void amrex_lo_averagecc (
        amrex_real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
        const amrex_real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
        const int *tlo, const int *thi, const int *nc
        );

    void amrex_lo_averageec (
        amrex_real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
        const amrex_real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
        const int *tlo, const int *thi, const int *nc,
        const int *axis
        );

    void amrex_lo_harmonic_averageec (
        amrex_real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
        const amrex_real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
        const int *tlo, const int *thi, const int *nc,
        const int *axis
        );
#ifdef __cplusplus
}
#endif

#endif /*_LO_F_H_*/
