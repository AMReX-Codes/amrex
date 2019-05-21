#ifndef WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_
#define WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_

#include <limits>
#include <WarpXParticleContainer.H>
#include <AMReX_REAL.H>

#ifndef WARPX_RZ

void
GetPosition(
    amrex::Real& x, amrex::Real& y, amrex::Real& z,
    const WarpXParticleContainer::ParticleType& p)
{
#if (AMREX_SPACEDIM==3)
    x = p.pos(0);
    y = p.pos(1);
    z = p.pos(2);
#else
    x = p.pos(0);
    z = p.pos(1);
#endif
}

# else // if WARPX_RZ is True

void
GetCartesianPositionFromCylindrical(
    amrex::Real& x, amrex::Real& y, amrex::Real& z,
    const WarpXParticleContainer::ParticleType& p, const amrex::Real theta)
{
    const amrex::Real r = p.pos(0);
    x = r*std::cos(theta);
    y = r*std::sin(theta);
    z = p.pos(1);
}

#endif // WARPX_RZ

#endif // WARPX_PARTICLES_PUSHER_GETANDSETPOSITION_H_
