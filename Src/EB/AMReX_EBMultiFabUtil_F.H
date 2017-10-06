#ifndef AMREX_EBMULTIFABUTIL_F_H_
#define AMREX_EBMULTIFABUTIL_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void amrex_eb_set_covered (const int* lo, const int* hi,
                               amrex_real* data, const int* dlo, const int* dhi,
                               const void* flag, const int* flo, const int* fhi,
                               const amrex_real* v, const int* ncomp);

    void amrex_eb_avgdown_sv (const int* lo, const int* hi,
                              const amrex_real* fine, const int* flo, const int* fhi,
                              amrex_real* crse, const int* clo, const int* chi,
                              const amrex_real* vol, const int* vlo, const int* vhi,
                              const amrex_real* vfol, const int* vflo, const int* vfhi,
                              const int* ratio, const int* ncomp);
    
#ifdef __cplusplus
}
#endif

#endif
