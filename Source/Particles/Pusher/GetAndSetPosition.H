#ifndef WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_
#define WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_

#include <limits>
#include <WarpXParticleContainer.H>
#include <AMReX_REAL.H>

#ifndef WARPX_RZ

/* \brief Extract the particle's coordinates from the ParticleType struct `p`,
 *        and stores them in the variables `x`, `y`, `z`. */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void GetPosition(
    amrex::Real& x, amrex::Real& y, amrex::Real& z,
    const WarpXParticleContainer::ParticleType& p)
{
#if (AMREX_SPACEDIM==3)
    x = p.pos(0);
    y = p.pos(1);
    z = p.pos(2);
#else
    x = p.pos(0);
    y = std::numeric_limits<amrex::Real>::quiet_NaN();
    z = p.pos(1);
#endif
}

/* \brief Set the particle's coordinates in the ParticleType struct `p`,
 *        from their values in the variables `x`, `y`, `z`. */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void SetPosition(
    WarpXParticleContainer::ParticleType& p,
    const amrex::Real x, const amrex::Real y, const amrex::Real z)
{
#if (AMREX_SPACEDIM==3)
    p.pos(0) = x;
    p.pos(1) = y;
    p.pos(2) = z;
#else
    p.pos(0) = x;
    p.pos(1) = z;
#endif
}

# else // if WARPX_RZ is True

/* \brief Extract the particle's coordinates from `theta` and the attributes
 *         of the ParticleType struct `p` (which contains the radius),
 *         and store them in the variables `x`, `y`, `z` */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void GetCartesianPositionFromCylindrical(
    amrex::Real& x, amrex::Real& y, amrex::Real& z,
    const WarpXParticleContainer::ParticleType& p, const amrex::Real theta)
{
    const amrex::Real r = p.pos(0);
    x = r*std::cos(theta);
    y = r*std::sin(theta);
    z = p.pos(1);
}

/* \brief Set the particle's cylindrical coordinates by setting `theta`
 *        and the attributes of the ParticleType struct `p` (which stores the radius),
 *        from the values of `x`, `y`, `z` */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void SetCylindricalPositionFromCartesian(
    WarpXParticleContainer::ParticleType& p, amrex::Real& theta,
    const amrex::Real x, const amrex::Real y, const amrex::Real z )
{
    theta = std::atan2(y, x);
    p.pos(0) = std::sqrt(x*x + y*y);
    p.pos(1) = z;
}

#endif // WARPX_RZ

#endif // WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_
