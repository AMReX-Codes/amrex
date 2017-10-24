#ifndef EBPARTICLES_F_H_
#define EBPARTICLES_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void amrex_move_particles(void* particles, const int* np,
                              const amrex::Real* dt,
                              const amrex::Real* problo,
                              const amrex::Real* problhi);

    void amrex_bounce_walls(const void* particles, const int* np,
                            const amrex_real* plo, const amrex_real* dx, 
                            const amrex::Real* dt,
                            const void* flag, const int* fglo, const int* fghi,
                            const amrex_real* bcent, const int* blo, const int* bhi,
                            const amrex_real* ax, const int* axlo, const int* axhi,
                            const amrex_real* ay, const int* aylo, const int* ayhi);

#ifdef __cplusplus
}
#endif

#endif /* _EBPARTICLES_F_H_ */
