#ifndef ELECTROMAGNETIC_PARTICLE_CONTAINER_H_
#define ELECTROMAGNETIC_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>

struct PIdx
{
    enum {
        x = 0,
        y, z, ux, uy, uz, w, Ex, Ey, Ez, Bx, By, Bz, ginv,
        nattribs
    };
};


class ElectromagneticParticleContainer
{

public:

    ElectromagneticParticleContainer (const amrex::Geometry            & a_geom,
                                      const amrex::DistributionMapping & a_dmap,
                                      const amrex::BoxArray            & a_ba,
                                      const int                         a_species_id,
                                      const amrex::Real                 a_charge,
                                      const amrex::Real                 a_mass);

    void InitParticles(const amrex::IntVect& a_num_particles_per_cell,
                       const amrex::Real     a_thermal_momentum_std,
                       const amrex::Real     a_thermal_momentum_mean,
                       const amrex::Real     a_density,
                       const amrex::RealBox& a_bounds,
                       const int             a_problem);

    void PushAndDeposeParticles(const amrex::MultiFab& Ex,
                                const amrex::MultiFab& Ey,
                                const amrex::MultiFab& Ez,
                                const amrex::MultiFab& Bx,
                                const amrex::MultiFab& By,
                                const amrex::MultiFab& Bz,
                                amrex::MultiFab& jx,
                                amrex::MultiFab& jy,
                                amrex::MultiFab& jz,
                                amrex::Real      dt);

    void PushParticleMomenta(const amrex::MultiFab& Ex,
                             const amrex::MultiFab& Ey,
                             const amrex::MultiFab& Ez,
                             const amrex::MultiFab& Bx,
                             const amrex::MultiFab& By,
                             const amrex::MultiFab& Bz,
                             amrex::Real      dt);

    void PushParticlePositions(amrex::Real dt);

    void EnforcePeriodicBCs();

protected:

    int m_species_id;

    amrex::Real m_charge;
    amrex::Real m_mass;
};

#endif
