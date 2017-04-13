#ifndef AMREX_COMPUTE_FORCE_F_H_
#define AMREX_COMPUTE_FORCE_F_H_

#include <AMReX_BLFort.H>

extern "C"
{

    void amrex_compute_forces(const amrex::Real*, int ns, int np,
                              amrex::Real* ghosts, int ng,
                              D_DECL(amrex::Real* ax,
                                     amrex::Real* ay,
                                     amrex::Real* az));

}

#endif /* _COMPUTE_FORCE_F_H_ */
