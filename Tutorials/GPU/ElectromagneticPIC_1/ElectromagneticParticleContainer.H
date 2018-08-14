#ifndef ELECTROMAGNETIC_PARTICLE_CONTAINER_H_
#define ELECTROMAGNETIC_PARTICLE_CONTAINER_H_

#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <AMReX_REAL.H>
#include <AMReX_MultiFab.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>

#include <thrust/device_vector.h>

#include <Particles.H>

// Kernel Functions
// ==========================
AMREX_CUDA_GLOBAL
void
PushAndDeposeParticles (const amrex::BaseFab<amrex::Real>    &elecX, const amrex::BaseFab<amrex::Real> &elecY, const amrex::BaseFab<amrex::Real> &elecZ,
                        const amrex::BaseFab<amrex::Real>     &magX, const amrex::BaseFab<amrex::Real>  &magY, const amrex::BaseFab<amrex::Real>  &magZ,
                        amrex::BaseFab<amrex::Real> &currDenX, amrex::BaseFab<amrex::Real> &currDenY, amrex::BaseFab<amrex::Real> &currDenZ,
                        amrex::GeometryData geomData, ParticlesData pData,
                        amrex::Real charge, amrex::Real mass, amrex::Real dt);

AMREX_CUDA_GLOBAL
void PushParticleMomenta (const amrex::BaseFab<amrex::Real> &elecX, const amrex::BaseFab<amrex::Real> &elecY, const amrex::BaseFab<amrex::Real> &elecZ,
                          const amrex::BaseFab<amrex::Real>  &magX, const amrex::BaseFab<amrex::Real>  &magY, const amrex::BaseFab<amrex::Real>  &magZ,
                          amrex::GeometryData geomData, ParticlesData pData, 
                          amrex::Real charge, amrex::Real mass, amrex::Real dt);

AMREX_CUDA_GLOBAL
void PushParticlePositions(ParticlesData pdata, amrex::Real dt);

AMREX_CUDA_GLOBAL
void EnforcePeriodicBCs (ParticlesData pData, amrex::GeometryData geomData);


// ==========================


class ElectromagneticParticleContainer
{
    amrex::BoxArray m_ba;
    amrex::Geometry m_geom;
    amrex::DistributionMapping m_dmap;

    std::map<int, Particles> m_particles;
    std::unique_ptr<amrex::iMultiFab> m_mask_ptr;

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

    void OK();

    void Redistribute();

protected:

    void RedistributeMPI(std::map<int, thrust::device_vector<char> >& not_ours);
    
    int m_species_id;
    
    amrex::Real m_charge;
    amrex::Real m_mass;

    size_t superparticle_size;
};

#endif
