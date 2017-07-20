#ifndef _NEIGHBORPARTICLES_H_
#define _NEIGHBORPARTICLES_H_

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
/// Note: For the neighbor particles, we don't communicate the integer components, only the
/// real data.
///
template <int NStructReal, int NStructInt, int NNeighborReal>
class NeighborParticleContainer
    : public ParticleContainer<NStructReal, NStructInt, 0, 0>
{
public:

    using ParticleType = typename ParticleContainer<NStructReal,
                                                    NStructInt, 0, 0>::ParticleType;
    using MyParIter = ParIter<NStructReal, NStructInt, 0, 0>;
    using PairIndex = std::pair<int, int>;
    using NeighborCommMap = std::map<NeighborCommTag, Array<char> >;
    using AoS = typename ParticleContainer<NStructReal, NStructInt, 0, 0>::AoS;

    NeighborParticleContainer(ParGDBBase* gdb, int ncells);

    NeighborParticleContainer(const Geometry            & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba,
                              int                         nneighbor);

    void InitLevelMask(int lev, 
                       const Geometry &geom, 
                       const DistributionMapping &dmap, 
                       const BoxArray &ba );

    ///
    /// This fills the neighbor buffers for each tile with the proper data
    ///
    void fillNeighbors(int lev);

    ///
    /// This updates the neighbors with their current particle data.
    ///
    void updateNeighbors(int lev);

    ///
    /// Each tile clears its neighbors, freeing the memory
    ///
    void clearNeighbors(int lev);

    ///
    /// Build a Neighbor List for each tile
    ///
    void buildNeighborList(int lev, bool sort=false);

    void buildNeighborListFort(int lev, bool sort=false);

    std::map<PairIndex, Array<char> > neighbors;
    std::map<PairIndex, Array<int> > neighbor_list;
    const size_t pdata_size = (NNeighborReal+BL_SPACEDIM)*
        sizeof(typename ParticleType::RealType);

protected:
    
    ///
    /// Apply periodic shift to particle position, so naive distance calculation
    /// between neighbors and regular particles will be correct.
    ///
    void applyPeriodicShift(int lev, ParticleType& p, const IntVect& neighbor_cell);

    ///
    /// Pack a particle's data into the proper neighbor buffer, or put it into
    /// the structure to be sent to the other processes
    ///
    void packNeighborParticle(int lev,
                              const IntVect& neighbor_cell,
                              const BaseFab<int>& mask,
                              const ParticleType& p,
                              NeighborCommMap& neighbors_to_comm);
    
    ///
    /// Perform the MPI communication neccesary to fill neighbor buffers
    ///
    void fillNeighborsMPI(NeighborCommMap& neighbors_to_comm);

    virtual bool check_pair(const ParticleType& p1, const ParticleType& p2) {
        return false;
    };

    int num_neighbor_cells;
    FabArray<BaseFab<int> > mask;

    void cacheNeighborParticle(int lev, const PairIndex& index, int particle_index,
                               const IntVect& neighbor_cell, const BaseFab<int>& mask);

    std::map<PairIndex, Array<int> > buffer_id_cache;
    std::map<PairIndex, Array<IntVect> > buffer_cell_cache;
};

#include "AMReX_NeighborParticlesI.H"

}

#endif // _NEIGHBORPARTICLES_H_
