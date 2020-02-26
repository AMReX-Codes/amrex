/* Copyright 2019 Andrew Myers, David Grote, Maxence Thevenet
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_RigidInjectedParticleContainer_H_
#define WARPX_RigidInjectedParticleContainer_H_

#include "PhysicalParticleContainer.H"
#include <AMReX_Vector.H>

/**
 * When injecting a particle beam (typically for a plasma wakefield
 * acceleration simulation), say propagating in the z direction, it can
 * be necessary to make particles propagate in a straight line up to a given
 * location z=z0. This is of particular importance when running in a boosted
 * frame, where the beam may evolve due to its space charge fields before
 * entering the plasma, causing the actual injected beam, and hence the whole
 * simulation result, to depend on the Lorentz factor of the boost.
 *
 * This feature is implemented in RigidInjectedParticleContainer: At each
 * iteration, for each particle, if z<z0 the particle moves in a straight line,
 * and if z>z0 the particle evolves as a regular PhysicalParticleContainer.
 *
 * Note: This option is also useful to build self-consistent space charge
 * fields for the particle beam.
 *
 * RigidInjectedParticleContainer derives from PhysicalParticleContainer.
 */
class RigidInjectedParticleContainer
    : public PhysicalParticleContainer
{
public:
    RigidInjectedParticleContainer (amrex::AmrCore* amr_core,
                                    int ispecies,
                                    const std::string& name);
    virtual ~RigidInjectedParticleContainer () {}

    virtual void InitData() override;

    virtual void RemapParticles();
    virtual void BoostandRemapParticles();

    virtual void Evolve (int lev,
                         const amrex::MultiFab& Ex,
                         const amrex::MultiFab& Ey,
                         const amrex::MultiFab& Ez,
                         const amrex::MultiFab& Bx,
                         const amrex::MultiFab& By,
                         const amrex::MultiFab& Bz,
                         amrex::MultiFab& jx,
                         amrex::MultiFab& jy,
                         amrex::MultiFab& jz,
                         amrex::MultiFab* cjx,
                         amrex::MultiFab* cjy,
                         amrex::MultiFab* cjz,
                         amrex::MultiFab* rho,
                         amrex::MultiFab* crho,
                         const amrex::MultiFab* cEx,
                         const amrex::MultiFab* cEy,
                         const amrex::MultiFab* cEz,
                         const amrex::MultiFab* cBx,
                         const amrex::MultiFab* cBy,
                         const amrex::MultiFab* cBz,
                         amrex::Real t,
                         amrex::Real dt,
                         DtType a_dt_type=DtType::Full) override;

    virtual void PushPX (WarpXParIter& pti, amrex::Real dt, DtType a_dt_type=DtType::Full) override;

    virtual void PushP (int lev, amrex::Real dt,
                        const amrex::MultiFab& Ex,
                        const amrex::MultiFab& Ey,
                        const amrex::MultiFab& Ez,
                        const amrex::MultiFab& Bx,
                        const amrex::MultiFab& By,
                        const amrex::MultiFab& Bz) override;

    virtual void ReadHeader (std::istream& is) override;

    virtual void WriteHeader (std::ostream& os) const override;

private:

    // User input quantities
    amrex::Real zinject_plane = 0.;
    bool projected = true; // When true, particle transverse positions are directly projected (without adjusment)
    bool focused = false; // When true, particle transverse positions are adjusted to account for distance between zinject and z=0
    bool rigid_advance = true; // When true, particles are advance with vzbar before injection

    amrex::Real vzbeam_ave_boosted;

    amrex::Vector<int> done_injecting;
    amrex::Vector<amrex::Real> zinject_plane_levels;

    // Temporary quantites
    amrex::Real zinject_plane_lev;
    amrex::Real zinject_plane_lev_previous;
    bool done_injecting_lev;

};

#endif
