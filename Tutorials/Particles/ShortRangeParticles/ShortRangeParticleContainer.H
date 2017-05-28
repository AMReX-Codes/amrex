#ifndef _SHORTRANGEPARTICLECONTAINER_H_
#define _SHORTRANGEPARTICLECONTAINER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"

namespace {

struct NeighborCommTag {

    NeighborCommTag(int pid, int gid, int tid)
        : proc_id(pid), grid_id(gid), tile_id(tid)
    {}

    int proc_id;
    int grid_id;
    int tile_id;
};

bool operator<(const NeighborCommTag& l, const NeighborCommTag& r) {
    return (l.proc_id < r.proc_id || 
           (l.proc_id == r.proc_id && l.grid_id < r.grid_id) ||
           (l.proc_id == r.proc_id && l.grid_id == r.grid_id && l.tile_id < r.tile_id ));
}

}

namespace amrex {

///
/// This is a container for particles that undergo short-range interations.
/// In addition to the normal ParticleContainer methods, each tile contains a "neighbor
/// buffer" that is filled with data corresponding to the particles within 1 neighbor cell of
/// the tile boundaries. This allows the N^2 search over each pair of particles to proceed
/// locally, instead of over the entire domain.
///
/// Note that neighbor particles are different than "ghost" particles, which are used
/// in AMR subcycling to keep track of coarse level particles that may move on to fine
/// levels during a fine level time step.
///
/// The particles in this container carry velocity and acceleration in addition to the 
/// usual components.
///
/// In 2D, the struct is: x y vx vy ax ay id cpu
/// In 3D: x y z vx vy vz ax ay az id cpu
///
/// Note: For the neighbor particles, we don't communicate the integer components, only the
/// real data.
///
class ShortRangeParticleContainer
    : public ParticleContainer<2*BL_SPACEDIM>
{
public:

    using MyParIter = amrex::ParIter<2*BL_SPACEDIM>;
    using PairIndex = std::pair<int, int>;
    using NeighborCommMap = std::map<NeighborCommTag, Array<char> >;
    
    ///
    /// This particle container fills a mask for quickly computing
    /// neighbor grids / tiles for a given particle
    ///
    ShortRangeParticleContainer(const Geometry            & geom, 
                                const DistributionMapping & dmap,
                                const BoxArray            & ba,
                                int                         nneighbor);

    ///
    /// Init one particle per cell with random velocities
    ///
    void InitParticles();

    ///
    /// This fills the neighbor buffers for each tile with the proper data
    ///
    void fillNeighbors();

    ///
    /// Each tile clears its neighbors, freeing the memory
    ///
    void clearNeighbors();

    ///
    /// Compute the short range forces on a tile's worth of particles.
    /// fillNeighbors must have already been called.
    ///
    void computeForces();

    ///
    /// Move the particles according to their forces, reflecting at domain boundaries
    ///
    void moveParticles(const Real dt);

    ///
    /// Save the particle data in an ASCII format
    ///
    void writeParticles(int n);

private:
    
    ///
    /// Apply periodic shift to particle position, so naive distance calculation
    /// between neighbors and regular particles will be correct.
    ///
    void applyPeriodicShift(ParticleType& p, const IntVect& neighbor_cell);

    ///
    /// Pack a particle's data into the proper neighbor buffer, or put it into
    /// the structure to be sent to the other processes
    ///
    void packNeighborParticle(const IntVect& neighbor_cell,
                              const BaseFab<int>& mask,
                              const ParticleType& p,
                              NeighborCommMap& neighbors_to_comm);
    
    ///
    /// Perform the MPI communication neccesary to fill neighbor buffers
    ///
    void fillNeighborsMPI(NeighborCommMap& neighbors_to_comm);

    const int lev = 0;
    int num_neighbor_cells;
    const size_t pdata_size = 3*BL_SPACEDIM*sizeof(ParticleType::RealType);
    FabArray<BaseFab<int> > mask;
    std::map<PairIndex, Array<char> > neighbors;
};

}

#endif // _SHORTRANGEPARTICLECONTAINER_H_
