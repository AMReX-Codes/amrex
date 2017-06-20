#ifndef AMREX_COMPUTE_FORCE_F_H_
#define AMREX_COMPUTE_FORCE_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void amrex_compute_forces(      void* particles, const int* np,
                              const void* ghosts, const int* nn);

    void amrex_move_particles(void* particles, const int* np,
                              const amrex::Real* dt,
                              const amrex::Real* problo,
                              const amrex::Real* problhi);

#ifdef __cplusplus
}
#endif

#endif /* _COMPUTE_FORCE_F_H_ */
