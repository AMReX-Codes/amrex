/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_BORIS_WITHRR_H_
#define WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_BORIS_WITHRR_H_

#include "UpdateMomentumBoris.H"

#include <AMReX_REAL.H>


/**
 * Push the particle's positions over one timestep,
 * given the value of its momenta `ux`, `uy`, `uz`.
 * Includes Radiation Reaction according to
 * https://doi.org/10.1088/1367-2630/12/12/123005
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdateMomentumBorisWithRadiationReaction(
    amrex::ParticleReal& ux, amrex::ParticleReal& uy, amrex::ParticleReal& uz,
    const amrex::ParticleReal Ex, const amrex::ParticleReal Ey, const amrex::ParticleReal Ez,
    const amrex::ParticleReal Bx, const amrex::ParticleReal By, const amrex::ParticleReal Bz,
    const amrex::Real q, const amrex::Real m, const amrex::Real dt )
{
    //RR algorithm needs to store old value of the normalized momentum
    const amrex::Real ux_old = ux;
    const amrex::Real uy_old = uy;
    const amrex::Real uz_old = uz;

    //Useful constant
    constexpr amrex::Real inv_c2 = 1./(PhysConst::c*PhysConst::c);

    //Call to regular Boris pusher
    UpdateMomentumBoris(
        ux, uy, uz,
        Ex, Ey, Ez,
        Bx, By, Bz,
        q, m, dt );

    //Estimation of the normalized momentum at intermediate (integer) time
    const amrex::Real ux_n = (ux+ux_old)*0.5;
    const amrex::Real uy_n = (uy+uy_old)*0.5;
    const amrex::Real uz_n = (uz+uz_old)*0.5;

    // Compute Lorentz factor (and inverse) at intermediate (integer) time
    const amrex::Real gamma_n = std::sqrt( 1. +
        (ux_n*ux_n + uy_n*uy_n + uz_n*uz_n)*inv_c2);
    const amrex::Real inv_gamma_n = 1.0/gamma_n;

    //Estimation of the velocity at intermediate (integer) time
    const amrex::Real vx_n = ux_n*inv_gamma_n;
    const amrex::Real vy_n = uy_n*inv_gamma_n;
    const amrex::Real vz_n = uz_n*inv_gamma_n;
    const amrex::Real bx_n = vx_n/PhysConst::c;
    const amrex::Real by_n = vy_n/PhysConst::c;
    const amrex::Real bz_n = vz_n/PhysConst::c;

    //Lorentz force over charge
    const amrex::Real flx_q = (Ex + vy_n*Bz - vz_n*By);
    const amrex::Real fly_q = (Ey + vz_n*Bx - vx_n*Bz);
    const amrex::Real flz_q = (Ez + vx_n*By - vy_n*Bx);
    const amrex::Real fl_q2 = flx_q*flx_q + fly_q*fly_q + flz_q*flz_q;

    //Calculation of auxiliary quantities
    const amrex::Real bdotE = (bx_n*Ex + by_n*Ey + bz_n*Ez);
    const amrex::Real bdotE2 = bdotE*bdotE;
    const amrex::Real coeff = gamma_n*gamma_n*(fl_q2-bdotE2);

    //Radiation reaction constant
    const amrex::Real RRcoeff = 2.0*PhysConst::r_e*q*q/
        (3.0*m*m*PhysConst::c*PhysConst::c);

    //Compute the components of the RR force
    const amrex::Real frx =
        RRcoeff*(PhysConst::c*(fly_q*Bz - flz_q*By) + bdotE*Ex - coeff*bx_n);
    const amrex::Real fry =
        RRcoeff*(PhysConst::c*(flz_q*Bx - flx_q*Bz) + bdotE*Ey - coeff*by_n);
    const amrex::Real frz =
        RRcoeff*(PhysConst::c*(flx_q*By - fly_q*Bx) + bdotE*Ez - coeff*bz_n);

    //Update momentum using the RR force
    ux += frx*dt;
    uy += fry*dt;
    uz += frz*dt;
}

#endif // WARPX_PARTICLES_PUSHER_UPDATEMOMENTUM_BORIS_WITHRR_H_
