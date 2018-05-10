#ifndef ELECTROMAGNETIC_PARTICLE_CONTAINER_H_
#define ELECTROMAGNETIC_PARTICLE_CONTAINER_H_

#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <AMReX_REAL.H>
#include <AMReX_MultiFab.H>
#include <AMReX_IntVect.H>

#include <Particles.H>

class ElectromagneticParticleContainer
{
    amrex::BoxArray m_ba;
    amrex::Geometry m_geom;
    amrex::DistributionMapping m_dmap;

    std::map<int, Particles> m_particles;
    std::unique_ptr<amrex::MultiFab> m_mask_ptr;

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
                       const amrex::Real     a_density);
    
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
    
    void PushParticlesOnly(amrex::Real dt);

    void EnforcePeriodicBCs();

    void OK();

    void Redistribute();

protected:
    
    int m_species_id;
    
    amrex::Real m_charge;
    amrex::Real m_mass;
};

#endif
