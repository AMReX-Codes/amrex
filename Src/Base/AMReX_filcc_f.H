#ifndef AMREX_FILLCC_F_H_
#define AMREX_FILLCC_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_fab_filcc (amrex_real* q, const int* qlo, const int* qhi, const int* nq,
                          const int* domlo, const int* domhi,
                          const amrex_real* dx, const amrex_real* xlo,
                          const int* bclo);


    void amrex_hoextraptocc (amrex_real* q, const int* qlo, const int* qhi,
                             const int* domlo, const int* domhi,
                             const amrex_real* dx, const amrex_real* xlo);

#ifdef __cplusplus
}
#endif

#endif
