#ifndef _CELL_SORTED_F_H_
#define _CELL_SORTED_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void move_particles(void* particles,
                        const int *np,
                        const int *lo,
                        const int *hi,
                        int** c_vectors,
                        int* sizes,
                        const int* clo,
                        const int* chi,
                        const amrex_real* plo,
                        const amrex_real* dx,
                        const amrex_real* dt);
    
#ifdef __cplusplus
};
#endif

#endif /*_EM_PIC_F_H_*/
