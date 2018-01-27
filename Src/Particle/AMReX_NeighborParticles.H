#ifndef _NEIGHBORPARTICLES_H_
#define _NEIGHBORPARTICLES_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"

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

    struct NeighborIndexMap {
        int dst_grid;
        int dst_tile;
        int dst_index;
        int src_grid;
        int src_tile;
        int src_index;
        int thread_num;

        NeighborIndexMap(int dgrid, int dtile, int dindex,
                         int sgrid, int stile, int sindex, int tnum)
            : dst_grid(dgrid), dst_tile(dtile), dst_index(dindex),
              src_grid(sgrid), src_tile(stile), src_index(sindex), thread_num(tnum)
        {}
    };
    
    struct NeighborCopyTag {
        int grid;
        int tile;
        int src_index;
        int dst_offset;
        int periodic_shift[3];

        bool operator<(const NeighborCopyTag& other) const {
            if (grid != other.grid) return grid < other.grid;
            if (tile != other.tile) return tile < other.tile;
            if (periodic_shift[0] != other.periodic_shift[0]) 
                return periodic_shift[0] < other.periodic_shift[0];
            if (periodic_shift[1] != other.periodic_shift[1]) 
                return periodic_shift[1] < other.periodic_shift[1];
            if (periodic_shift[2] != other.periodic_shift[2]) 
                return periodic_shift[2] < other.periodic_shift[2];
            return false;
        }
        
        bool operator==(const NeighborCopyTag& other) const {
            return (grid == other.grid) and (tile == other.tile)
                and (periodic_shift[0] == other.periodic_shift[0])
                and (periodic_shift[1] == other.periodic_shift[1])
                and (periodic_shift[2] == other.periodic_shift[2]);
        }
    };

    struct NeighborCommTag {
        
        NeighborCommTag(int pid, int gid, int tid)
            : proc_id(pid), grid_id(gid), tile_id(tid)
            {}
        
        int proc_id;
        int grid_id;
        int tile_id;
        
        bool operator<(const NeighborCommTag& other) const {
            return (proc_id < other.proc_id || 
                    (proc_id == other.proc_id && 
                     grid_id < other.grid_id) ||
                    (proc_id == other.proc_id && 
                     grid_id == other.grid_id && 
                     tile_id < other.tile_id ));
        }

        bool operator==(const NeighborCommTag& other) const {
            return ( (proc_id == other.proc_id) and 
                     (grid_id == other.grid_id) and
                     (tile_id == other.tile_id) );
        }
    };    
        
public:

    using ParticleType = typename ParticleContainer<NStructReal,
                                                    NStructInt, 0, 0>::ParticleType;
    using MyParIter = ParIter<NStructReal, NStructInt, 0, 0>;
    using PairIndex = std::pair<int, int>;
    using NeighborCommMap = std::map<NeighborCommTag, Vector<char> >;
    using AoS = typename ParticleContainer<NStructReal, NStructInt, 0, 0>::AoS;

    NeighborParticleContainer(ParGDBBase* gdb, int ncells);

    NeighborParticleContainer(const Geometry            & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba,
                              int                         nneighbor);

    ///
    /// This resets the particle container to use the given BoxArray
    /// and DistributionMapping
    ///
    void Regrid(const DistributionMapping &dmap, const BoxArray &ba );

    ///
    /// This builds the internal data structure used for looking up neighbors
    ///
    void BuildLevelMask(int lev);

    ///
    /// This fills the neighbor buffers for each tile with the proper data
    ///
    void fillNeighbors(int lev);

    void cacheNeighborInfo(int lev);

    ///
    /// This updates the neighbors with their current particle data.
    ///
    void updateNeighbors(int lev, bool reuse_rcv_counts=true);

    ///
    /// Each tile clears its neighbors, freeing the memory
    ///
    void clearNeighbors(int lev);

    ///
    /// Build a Neighbor List for each tile
    ///
    void buildNeighborList(int lev, bool sort=false);

    void buildNeighborListFort(int lev, bool sort=false);

    void setRealCommComp(int i, bool value);
    void setIntCommComp(int i, bool value);

    std::map<PairIndex, Vector<char> > neighbors;
    std::map<PairIndex, Vector<int>  > neighbor_list;
    const size_t pdata_size = sizeof(ParticleType);
    
protected:
    
    void initializeCommComps();
    
    void calcCommSize();
    
    ///
    /// Perform the MPI communication neccesary to fill neighbor buffers
    ///
    void fillNeighborsMPI(bool reuse_rcv_counts);

    ///
    /// Perform handshake to figure out how many bytes each proc should receive
    ///
    void getRcvCountsMPI();

    virtual bool check_pair(const ParticleType& p1, const ParticleType& p2) {
        return false;
    };

    size_t cdata_size;
    int num_neighbor_cells;
    amrex::Vector<NeighborCommTag> local_neighbors;
    std::unique_ptr<FabArray<BaseFab<int> > > mask_ptr;

    std::map<PairIndex, Vector<Vector<NeighborCopyTag> > > buffer_tag_cache;
    std::map<PairIndex, int> local_neighbor_sizes;

    // each proc knows how many sends it will do, and how many bytes it will rcv 
    // from each other proc.
    amrex::Vector<int> neighbor_procs;
    amrex::Vector<long> rcvs;
    long num_snds;
    std::map<int, Vector<char> > send_data;

    std::array<bool, BL_SPACEDIM + NStructReal> rc; 
    std::array<bool, 2 + NStructInt>  ic;
};

#include "AMReX_NeighborParticlesI.H"

}

#endif // _NEIGHBORPARTICLES_H_
