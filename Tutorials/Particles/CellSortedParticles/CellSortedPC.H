#ifndef CELL_SORTED_PARTICLE_CONTAINER_H_
#define CELL_SORTED_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>

struct RealData
{
    enum {
        vx = 0,
        vy,
        vz,
        ncomps
    };
};

struct IntData
{
    enum {
        sorted = 0,
        ncomps
    };
};

class CellSortedParticleContainer
    : public amrex::ParticleContainer<RealData::ncomps, IntData::ncomps>
{
public:

    using MyParIter = amrex::ParIter<RealData::ncomps, IntData::ncomps, 0, 0>;
    
    CellSortedParticleContainer (const amrex::Geometry            & a_geom,
                                 const amrex::DistributionMapping & a_dmap,
                                 const amrex::BoxArray            & a_ba);

    void InitParticles(const amrex::IntVect& a_num_particles_per_cell);

    void InitCellVectors();
    
protected:
    
    std::map<int, amrex::BaseFab<std::unique_ptr<std::vector<int> > > > m_cell_vectors;
    std::map<int, amrex::BaseFab<int> > m_vector_size;
    std::map<int, amrex::BaseFab<int*> > m_vector_ptrs;
};

#endif
