#ifndef TEST_PARTICLE_CONTAINER_H_
#define TEST_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <AMReX_REAL.H>
#include <AMReX_MultiFab.H>
#include <AMReX_IntVect.H>

struct PIdx
{
    enum { // These are "permanent" attributes that move around with the particles
        vx = 0,
        vy, vz, w,
        nattribs
    };
};

class TestParticleContainer
    : public amrex::ParticleContainer<0, 0, PIdx::nattribs>
{
    
public:

    using TestParIter = amrex::ParIter<0, 0, PIdx::nattribs>;
    
    TestParticleContainer (const amrex::Geometry            & a_geom,
                           const amrex::DistributionMapping & a_dmap,
                           const amrex::BoxArray            & a_ba);
    
    void InitParticles(const amrex::IntVect& a_num_particles_per_cell,
                       const amrex::Real     a_thermal_momentum_std,
                       const amrex::Real     a_thermal_momentum_mean,
                       const amrex::Real     a_density,
                       const amrex::RealBox& a_bounds);
    
    void MoveParticles();
};

#endif
