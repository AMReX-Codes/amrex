#ifndef WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_BORIS_WITHRR_H_
#define WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_BORIS_WITHRR_H_

#include <AMReX_REAL.H>
#include <UpdateMomentumBoris.H>

/* \brief Push the particle's positions over one timestep,
 *    given the value of its momenta `ux`, `uy`, `uz`.
 *   Includes RR according to https://doi.org/10.1088/1367-2630/12/12/123005 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdateMomentumBorisWithRadiationReaction(
    amrex::Real& ux, amrex::Real& uy, amrex::Real& uz,
    const amrex::Real Ex, const amrex::Real Ey, const amrex::Real Ez,
    const amrex::Real Bx, const amrex::Real By, const amrex::Real Bz,
    const amrex::Real q, const amrex::Real m, const amrex::Real dt )
{
    //RR algorithm needs to store old value of momentum
    const amrex::Real ux_old = ux;
    const amrex::Real uy_old = uy;
    const amrex::Real uz_old = uz;

    //Useful constants
    const amrex::Real inv_dt = 1.0/dt;
    constexpr amrex::Real inv_c2 = 1./(PhysConst::c*PhysConst::c);

    //Call to regular Boris pusher
    UpdateMomentumBoris(
        ux, uy, uz,
        Ex, Ey, Ez,
        Bx, By, Bz,
        q, m, dt );

    //Estimation of the momentum at intermediate (integer) time
    const amrex::Real ux_n = (ux-ux_old)*0.5;
    const amrex::Real uy_n = (uy-uy_old)*0.5;
    const amrex::Real uz_n = (uz-uz_old)*0.5;

    // Compute Lorentz factor (and inverse) at intermediate (integer) time
    const amrex::Real gamma_n = std::sqrt( 1. +
        (ux_n*ux_n + uy_n*uy_n + uz_n*uz_n)*inv_c2);
    const amrex::Real inv_gamma_n = 1.0/gamma_n;

    //Estimation of the velocity at intermediate (integer) time
    const amrex::Real vx_n = ux_n*inv_gamma;
    const amrex::Real vy_n = uy_n*inv_gamma;
    const amrex::Real vz_n = uz_n*inv_gamma;

    //Lorentz force
    const amrex::Real flx = q*(Ex + vy_n*Bz - vz_n*By);
    const amrex::Real fly = q*(Ey + vz_n*Bx - vx_n*Bz);
    const amrex::Real flz = q*(Ez + vx_n*By - vy_n*Bx);
    const amrex::Real fl2 = flx*flx + fly*fly + flz*flz;

    //TO DO !
    const amrex::Real RRcoeff = 1.0; //WRONG! SHOULD PUT THE RIGHT VALUE
    const amrex::Real frx = 0.0; //WRONG! SHOULD PUT THE RIGHT VALUE
    const amrex::Real fry = 0.0; //WRONG! SHOULD PUT THE RIGHT VALUE
    const amrex::Real frz = 0.0; //WRONG! SHOULD PUT THE RIGHT VALUE
    //TO DO !

    //Update momentum using the RR force
    ux -= frx*dt;
    uy -= fry*dt;
    uz -= frz*dt;
}

#endif // WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_BORIS_WITHRR_H_
