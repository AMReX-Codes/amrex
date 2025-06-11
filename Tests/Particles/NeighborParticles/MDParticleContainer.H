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
    : public amrex::NeighborParticleContainer<PIdx::ncomps, 1, 1, 1>
{

public:

    MDParticleContainer (const amrex::Geometry            & a_geom,
                         const amrex::DistributionMapping & a_dmap,
                         const amrex::BoxArray            & a_ba,
                         int                                a_numcells)
        : NeighborParticleContainer<PIdx::ncomps, 1, 1, 1>(a_geom, a_dmap, a_ba, a_numcells)
    {
        int num_runtime_real = 1;
        int num_runtime_int  = 1;
        for (int i = 0; i < num_runtime_real; ++i)
        {
            AddRealComp(true);
        }
        for (int i = 0; i < num_runtime_int; ++i)
        {
            AddIntComp(true);
        }
    }

    void InitParticles (const amrex::IntVect& a_num_particles_per_cell,
                        amrex::Real           a_thermal_momentum_std,
                        amrex::Real           a_thermal_momentum_mean);

    void writeParticles (int n);

    void reset_test_id ();

    void checkNeighborParticles ();

    void checkNeighborList ();

    std::pair<amrex::Real, amrex::Real>  minAndMaxDistance ();

    void moveParticles (amrex::ParticleReal dx);
};

#endif
