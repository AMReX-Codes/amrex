#ifndef MD_PARTICLE_CONTAINER_H_
#define MD_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

struct PIdx
{
    enum {
        vx = 0,
        vy, vz, ax, ay, az,
        ncomps
    };
};

class MDParticleContainer
    : public amrex::NeighborParticleContainer<PIdx::ncomps, 0>
{

public:

    MDParticleContainer (const amrex::Geometry            & a_geom,
                         const amrex::DistributionMapping & a_dmap,
                         const amrex::BoxArray            & a_ba,
                         int                                a_numcells)
        : NeighborParticleContainer<PIdx::ncomps, 0>(a_geom, a_dmap, a_ba, a_numcells)
    {}

    void InitParticles (const amrex::IntVect& a_num_particles_per_cell,
                        const amrex::Real     a_thermal_momentum_std,
                        const amrex::Real     a_thermal_momentum_mean);

    void computeForces ();

    void moveParticles (const amrex::Real& dt);

    void writeParticles (const int n);    
};

#endif
