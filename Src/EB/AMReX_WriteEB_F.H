#ifndef AMREX_EBWRITE_F_H_
#define AMREX_EBWRITE_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void amrex_eb_to_polygon(
            const amrex::Real* dx,
            const int* slo, const int* shi,
            const void* flag,        const int* fglo, const int* fghi,
            const amrex_real* bcent, const int* blo,  const int* bhi,
            const amrex_real* ax,    const int* axlo, const int* axhi,
            const amrex_real* ay,    const int* aylo, const int* ayhi,
            const amrex_real* az,    const int* azlo, const int* azhi
        );

    void amrex_write_eb_vtp(int* myID);
    void amrex_write_pvtp(int* nProcs);

    void amrex_eb_grid_coverage(
            int* myID, const amrex::Real* dx,
            const int* slo, const int* shi,
            const void* flag, const int* fglo, const int* fghi
        );

#ifdef __cplusplus
}
#endif

#endif
