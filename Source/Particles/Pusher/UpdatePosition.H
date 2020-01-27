#ifndef WARPX_PARTICLES_PUSHER_UPDATEPOSITION_H_
#define WARPX_PARTICLES_PUSHER_UPDATEPOSITION_H_

#include <AMReX_FArrayBox.H>
#include <WarpXConst.H>
#include <AMReX_REAL.H>

/** \brief Push the particle's positions over one timestep,
 *    given the value of its momenta `ux`, `uy`, `uz` */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdatePosition(WarpXParticleContainer::ParticleType& p,
                    const amrex::ParticleReal ux, const amrex::ParticleReal uy, const amrex::ParticleReal uz,
                    const amrex::Real dt )
{

    constexpr amrex::Real inv_c2 = 1./(PhysConst::c*PhysConst::c);

    // Compute inverse Lorentz factor
    const amrex::Real inv_gamma = 1./std::sqrt(1. + (ux*ux + uy*uy + uz*uz)*inv_c2);
    // Update positions over one time step
    p.pos(0) += ux * inv_gamma * dt;
#if (AMREX_SPACEDIM == 3) || (defined WARPX_DIM_RZ) // RZ pushes particles in 3D
    p.pos(1) += uy * inv_gamma * dt;
#endif
    p.pos(2) += uz * inv_gamma * dt;
}

#endif // WARPX_PARTICLES_PUSHER_UPDATEPOSITION_H_
