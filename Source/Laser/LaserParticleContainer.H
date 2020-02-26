/* Copyright 2019 Andrew Myers, Axel Huebl, David Grote
 * Luca Fedeli, Maxence Thevenet, Remi Lehe
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_LaserParticleContainer_H_
#define WARPX_LaserParticleContainer_H_

#include "LaserProfiles.H"

#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXConst.H"
#include "Parser/WarpXParser.H"

#include <memory>
#include <limits>

/**
 * The main method to inject a laser pulse in WarpX is to use an artificial
 * antenna: particles evenly distributed in a given plane (one particle per
 * cell) move at each iteration and deposit a current J onto the grid, which
 * in turns creates an electromagnetic field on the grid. The particles'
 * displacements are prescribed to create the field requested by the user.
 *
 * These artificial particles are contained in the LaserParticleContainer.
 * LaserParticleContainer derives directly from WarpXParticleContainer. It
 * requires a DepositCurrent function, but no FieldGather function.
 */
class LaserParticleContainer
    : public WarpXParticleContainer
{
public:
    LaserParticleContainer (amrex::AmrCore* amr_core, int ispecies, const std::string& name);
    virtual ~LaserParticleContainer () {}

    virtual void InitData () final;

#ifdef WARPX_DO_ELECTROSTATIC
    virtual void EvolveES (const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                           amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                           amrex::Real t, amrex::Real dt) { BL_ASSERT(false); }
#endif // WARPX_DO_ELECTROSTATIC

    virtual void Evolve (int lev,
                         const amrex::MultiFab&, const amrex::MultiFab&, const amrex::MultiFab&,
                         const amrex::MultiFab&, const amrex::MultiFab&, const amrex::MultiFab&,
                         amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz,
                         amrex::MultiFab*, amrex::MultiFab*, amrex::MultiFab*,
                         amrex::MultiFab* rho, amrex::MultiFab* crho,
                         const amrex::MultiFab*, const amrex::MultiFab*, const amrex::MultiFab*,
                         const amrex::MultiFab*, const amrex::MultiFab*, const amrex::MultiFab*,
                         amrex::Real t, amrex::Real dt, DtType a_dt_type=DtType::Full) final;

    virtual void PushP (int lev, amrex::Real dt,
                        const amrex::MultiFab& ,
                        const amrex::MultiFab& ,
                        const amrex::MultiFab& ,
                        const amrex::MultiFab& ,
                        const amrex::MultiFab& ,
                        const amrex::MultiFab& ) final;

    virtual void PostRestart () final;

    void calculate_laser_plane_coordinates (const WarpXParIter& pti, const int np,
                                            amrex::Real * AMREX_RESTRICT const pplane_Xp,
                                            amrex::Real * AMREX_RESTRICT const pplane_Yp);

    void update_laser_particle (WarpXParIter& pti, const int np, amrex::ParticleReal * AMREX_RESTRICT const puxp,
                                amrex::ParticleReal * AMREX_RESTRICT const puyp,
                                amrex::ParticleReal * AMREX_RESTRICT const puzp,
                                amrex::ParticleReal const * AMREX_RESTRICT const pwp,
                                amrex::Real const * AMREX_RESTRICT const amplitude,
                                const amrex::Real dt);

protected:

    std::string laser_name;

private:
    // runtime paramters
    amrex::Vector<amrex::Real> position; //! Coordinates of one of the point of the antenna
    amrex::Vector<amrex::Real> nvec; //! Normal of the plane of the antenna
    amrex::Vector<amrex::Real> p_X;// ! Polarization

    long                      pusher_algo = -1;
    amrex::Real               e_max       = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real               wavelength  = std::numeric_limits<amrex::Real>::quiet_NaN();

    amrex::Real               Z0_lab = 0; // Position of the antenna in the lab frame

    long min_particles_per_mode = 4;

    // computed using runtime parameters
    amrex::Vector<amrex::Real> p_Y;
    amrex::Vector<amrex::Real> u_X;
    amrex::Vector<amrex::Real> u_Y;
    amrex::Real weight   = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real mobility = std::numeric_limits<amrex::Real>::quiet_NaN();


    // laser particle domain
    amrex::RealBox laser_injection_box;
    // Theoretical position of the antenna. Used if do_continuous_injection=1.
    // Track the position of the antenna until it enters the simulation domain.
    amrex::Vector<amrex::Real> updated_position;

    void ComputeSpacing (int lev, amrex::Real& Sx, amrex::Real& Sy) const;
    void ComputeWeightMobility (amrex::Real Sx, amrex::Real Sy);
    void InitData (int lev);
    // Inject the laser antenna during the simulation, if it started
    // outside of the simulation domain and enters it.
    void ContinuousInjection(const amrex::RealBox& injection_box) override;
    // Update position of the antenna
    void UpdateContinuousInjectionPosition(amrex::Real dt) override;

    // Unique (smart) pointer to the laser profile
    std::unique_ptr<WarpXLaserProfiles::ILaserProfile> m_up_laser_profile;

};

#endif
