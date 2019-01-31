#ifndef AMREX_EBAMRUTIL_F_H_
#define AMREX_EBAMRUTIL_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_tag_cutcells (const int* lo, const int* hi,
                             char* tag, const int* tlo, const int* thi,
                             const void* flag, const int* flo, const int* fhi,
                             char tagval, char clearval);

    void amrex_tag_volfrac  (const int* lo, const int* hi,
                             char* tag, const int* tlo, const int* thi,
                             const void* vfrac, const int* vlo, const int* vhi,
                             char tagval, char clearval, amrex::Real tol);


#ifdef __cplusplus
}
#endif

#endif
