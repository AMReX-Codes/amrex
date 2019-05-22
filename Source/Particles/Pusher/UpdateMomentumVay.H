#ifndef WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_VAY_H_
#define WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_VAY_H_

#include <AMReX_FArrayBox.H>
#include <WarpXConst.H>
#include <AMReX_REAL.H>

/* \brief Push the particle's positions over one timestep,
 *    given the value of its momenta `ux`, `uy`, `uz` */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdateMomentumVay(
    amrex::Real& ux, amrex::Real& uy, amrex::Real& uz,
    const amrex::Real Ex, const amrex::Real Ey, const amrex::Real Ez,
    const amrex::Real Bx, const amrex::Real By, const amrex::Real Bz,
    const amrex::Real q, const amrex::Real m, const amrex::Real dt )
{
    // Constants
    const amrex::Real econst = q*dt/m;
    const amrex::Real bconst = 0.5*q*dt/m;
    constexpr amrex::Real inv_c2 = 1./(PhysConst::c*PhysConst::c);
    // Compute temporary gamma factor
    const amrex::Real inv_gamma = 1./std::sqrt(1. + (ux*ux + uy*uy + uz*uz)*inv_c2);
    // Compute temporary variables
    const amrex::Real taux = econst*Bx;
    const amrex::Real tauy = econst*By;
    const amrex::Real tauz = econst*Bz;
    const amrex::Real tausq = taux*taux + tauy*tauy + tauz*tauz;
    const amrex::Real ux_p = ux + econst*Ex + (uy*tauz - uz*tauy)*inv_gamma;
    const amrex::Real uy_p = uy + econst*Ey + (uz*taux - ux*tauz)*inv_gamma;
    const amrex::Real uz_p = uz + econst*Ez + (ux*tauy - uy*taux)*inv_gamma;
    // - Update momentum
    ux += uy_p*sz - uz_p*sy;
    uy += uz_p*sx - ux_p*sz;
    uz += ux_p*sy - uy_p*sx;
    // Second half-push for E
    ux += econst*Ex;
    uy += econst*Ey;
    uz += econst*Ez;
}

#endif // WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_VAY_H_
