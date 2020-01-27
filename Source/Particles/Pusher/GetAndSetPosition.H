#ifndef WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_
#define WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_

#include <limits>
#include <WarpXParticleContainer.H>
#include <AMReX_REAL.H>

#ifdef WARPX_DIM_RZ

/** \brief Extract the particle's coordinates from `theta` and the attributes
 *         of the ParticleType struct `p` (which contains the radius),
 *         and store them in the variables `x`, `y`, `z` */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void GetCartesianPositionFromCylindrical(
    amrex::ParticleReal& x, amrex::ParticleReal& y, amrex::ParticleReal& z,
    const WarpXParticleContainer::ParticleType& p, const amrex::ParticleReal theta)
{
    const amrex::ParticleReal r = p.pos(0);
    x = r*std::cos(theta);
    y = r*std::sin(theta);
    z = p.pos(1);
}

/** \brief Set the particle's cylindrical coordinates by setting `theta`
 *        and the attributes of the ParticleType struct `p` (which stores the radius),
 *        from the values of `x`, `y`, `z` */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void SetCylindricalPositionFromCartesian(
    WarpXParticleContainer::ParticleType& p, amrex::ParticleReal& theta,
    const amrex::ParticleReal x, const amrex::ParticleReal y, const amrex::ParticleReal z )
{
    theta = std::atan2(y, x);
    p.pos(0) = std::sqrt(x*x + y*y);
    p.pos(1) = z;
}

#endif // WARPX_DIM_RZ

#endif // WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_
