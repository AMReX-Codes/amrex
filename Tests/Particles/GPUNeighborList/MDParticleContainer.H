#ifndef MD_PARTICLE_CONTAINER_H_
#define MD_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>

struct PIdx
{
    enum {
        vx = 0,
        vy, vz, ax, ay, az,
        ncomps
    };
};

class MDParticleContainer
    : public amrex::ParticleContainer<PIdx::ncomps>
{

public:

    MDParticleContainer (const amrex::Geometry            & a_geom,
                         const amrex::DistributionMapping & a_dmap,
                         const amrex::BoxArray            & a_ba);

    void InitParticles(const amrex::IntVect& a_num_particles_per_cell,
                       const amrex::Real     a_thermal_momentum_std,
                       const amrex::Real     a_thermal_momentum_mean);

    void RedistributeLocal()
    {
        const int lev_min = 0;
        const int lev_max = 0;
        const int nGrow = 0;
        const int local = 1;
        Redistribute(lev_min, lev_max, nGrow, local);
    }

    void BuildNeighborList();

    void printNeighborList();

    void computeForces();

    void moveParticles(const amrex::Real& dt);

protected:
    
    // This is the neighbor list data structure
    std::map<std::pair<int, int>, amrex::Gpu::DeviceVector<unsigned int> > m_nbor_offsets;
    std::map<std::pair<int, int>, amrex::Gpu::DeviceVector<unsigned int> > m_nbor_list;
};

#endif
