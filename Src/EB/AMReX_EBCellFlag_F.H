#ifndef AMREX_EBCELLFLAG_F_H_
#define AMREX_EBCELLFLAG_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_ebcellflag_count(const int* lo, const int* hi,
                                const void* flag, const int* flo, const int* fhi,
                                int* nregular, int* nsingle, int* nmulti, int* ncovered);

#ifdef __cplusplus
}
#endif

#endif
