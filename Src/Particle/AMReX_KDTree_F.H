#ifndef AMREX_KDTREE_F_H_
#define AMREX_KDTREE_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void amrex_compute_best_partition(const amrex_real*, const int*, const int*,
                                      const int*, const int*, amrex_real, int,
                                      amrex_real*, amrex_real*, int*);

    void amrex_compute_cost(const amrex_real* pcounts, amrex_real* cost, const int* lo,
                            const int* hi, const amrex_real cell_weight);

#ifdef __cplusplus
}
#endif

#endif /* _AMREX_KDTREE_F_H_ */
