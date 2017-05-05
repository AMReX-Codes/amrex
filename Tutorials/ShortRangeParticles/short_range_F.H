#ifndef AMREX_COMPUTE_FORCE_F_H_
#define AMREX_COMPUTE_FORCE_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void amrex_compute_forces(const amrex_real*, int ns, int np,
                              const amrex_real* ghosts, int ng);

    void amrex_move_particles(const amrex_real*, int ns, int np,
                              const amrex_real* dt,
                              const amrex_real* problo,
                              const amrex_real* problhi);
#ifdef __cplusplus
}
#endif

#endif /* _COMPUTE_FORCE_F_H_ */
