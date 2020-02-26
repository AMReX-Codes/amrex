/* Copyright 2019 David Grote, Luca Fedeli, Maxence Thevenet
 * Remi Lehe, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_PUSHER_UPDATEPOSITIONPHOTON_H_
#define WARPX_PARTICLES_PUSHER_UPDATEPOSITIONPHOTON_H_

#include "Utils/WarpXConst.H"

#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>


/**
 * \brief Push the position of a photon particle over one timestep,
 *    given the value of its momenta `ux`, `uy`, `uz`
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdatePositionPhoton(
    amrex::ParticleReal& x, amrex::ParticleReal& y, amrex::ParticleReal& z,
    const amrex::ParticleReal ux, const amrex::ParticleReal uy, const amrex::ParticleReal uz,
    const amrex::Real dt )
{
    // Compute speed of light over inverse of momentum modulus
    const amrex::Real c_over_umod = PhysConst::c/std::sqrt(ux*ux + uy*uy + uz*uz);

    // Update positions over one time step
    x += ux * c_over_umod * dt;
#if (defined WARPX_DIM_3D) || (defined WARPX_DIM_RZ) // RZ pushes particles in 3D
    y += uy * c_over_umod * dt;
#endif
    z += uz * c_over_umod * dt;
}

#endif // WARPX_PARTICLES_PUSHER_UPDATEPOSITION_H_
