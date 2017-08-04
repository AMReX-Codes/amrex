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
    
#ifdef __cplusplus
}
#endif

#endif
