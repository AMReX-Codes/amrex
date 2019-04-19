#ifndef POISSON_F_H_
#define POISSON_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    
    // these are interfaces to fortran subroutines

    void init_data (const int* lo, const int* hi,
                    amrex_real* u, const int* ulo, const int* uhi);


#ifdef __cplusplus
}
#endif

#endif
