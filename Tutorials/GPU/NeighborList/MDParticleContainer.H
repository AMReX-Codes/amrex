#ifndef MD_PARTICLE_CONTAINER_H_
#define MD_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>
#include <AMReX_CudaContainers.H>

#include <AMReX_NeighborList.H>

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

    void InitParticles (const amrex::IntVect& a_num_particles_per_cell,
                        const amrex::Real     a_thermal_momentum_std,
                        const amrex::Real     a_thermal_momentum_mean);

    void RedistributeLocal ();

    void BuildNeighborList ();

    void printNeighborList ();

    void computeForces ();

    void moveParticles (const amrex::Real& dt);

    void sortParticlesByNeighborDest ();

    void fillNeighbors ();

    void clearNeighbors ();

    void writeParticles (const int n);

private:

    void fillNeighborsMPIGPU(std::map<int, SendBuffer>& not_ours);
 
    void buildNeighborMask ();

    // These are used to keep track of which particles need to be ghosted to which grids
    bool m_neighbor_mask_initialized = false;
    std::unique_ptr<amrex::iMultiFab> m_neighbor_mask_ptr;
    std::map<int, std::vector<std::vector<int> > > m_grid_map;
    std::map<int, amrex::Gpu::HostVector<int> > m_start;
    std::map<int, amrex::Gpu::HostVector<int> > m_stop;
    
    std::map<std::pair<int, int>, amrex::NeighborList<ParticleType> > m_neighbor_list;
};

#endif
