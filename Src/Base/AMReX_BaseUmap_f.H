#ifndef AMREX_BASEUMAP_F_H_
#define AMREX_BASEUMAP_F_H_

#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C"
{
#endif
 
    amrex_real amrex_fort_umap_norm (const int* lo, const int* hi,
                              const amrex_real* src, const int* src_sz,
                              const int* kt, const int* ktlo, const int* kthi, 
                              const int* max_mv, const int* ncomp,
                              const int* p);


#ifdef __cplusplus
}
#endif

#endif
