#ifndef AMREX_MG_F_H_
#define AMREX_MG_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void amrex_mg_average (
        amrex_real* crse,       AMREX_ARLIM_P(crse_lo), AMREX_ARLIM_P(crse_hi),
        const amrex_real* fine, AMREX_ARLIM_P(fine_lo), AMREX_ARLIM_P(fine_hi),
        const int *tlo, const int *thi,
        const int *nc);

    void amrex_mg_interp (
        amrex_real* fine,       AMREX_ARLIM_P(fine_lo), AMREX_ARLIM_P(fine_hi),
        const amrex_real* crse, AMREX_ARLIM_P(crse_lo), AMREX_ARLIM_P(crse_hi),
        const int *tlo, const int *thi,
        const int *nc);
#ifdef __cplusplus
}
#endif

#endif /*_MG_F_H_*/
