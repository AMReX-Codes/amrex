#ifndef AMREX_NEIGHBORPARTICLES_H_
#define AMREX_NEIGHBORPARTICLES_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>

#ifdef AMREX_USE_CUDA
#include <AMReX_NeighborList.H>
#endif

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
template <int NStructReal, int NStructInt>
class NeighborParticleContainer
    : public ParticleContainer<NStructReal, NStructInt, 0, 0>
{

    struct MaskComps
    {
        enum {
            grid = 0,
            tile,
            level
        };
    };

    struct NeighborIndexMap {
        int dst_level;
        int dst_grid;
        int dst_tile;
        int dst_index;
        int src_level;
        int src_grid;
        int src_tile;
        int src_index;
        int thread_num;

        NeighborIndexMap(int dlevel, int dgrid, int dtile, int dindex,
                         int slevel, int sgrid, int stile, int sindex, int tnum)
            : dst_level(dlevel), dst_grid(dgrid), dst_tile(dtile), dst_index(dindex),
              src_level(slevel), src_grid(sgrid), src_tile(stile), src_index(sindex),
              thread_num(tnum)
        {}

        friend std::ostream& operator<< (std::ostream& os, const NeighborIndexMap& nim)
        {
            os << nim.dst_level << " " << nim.dst_grid << " " << nim.dst_tile << " " << nim.dst_index
               << nim.src_level << " " << nim.src_grid << " " << nim.src_tile << " " << nim.src_index
               << nim.thread_num;
            if (!os.good())
                amrex::Error("operator<<(ostream&,Particle<NReal, NInt>&) failed");
            return os;
        }

    };

    struct NeighborCopyTag {
        int level;
        int grid;
        int tile;
        int src_index;
        int dst_index;
        IntVect periodic_shift;

        NeighborCopyTag() {}

        NeighborCopyTag(int a_level, int a_grid, int a_tile) :
            level(a_level), grid(a_grid), tile(a_tile), src_index(0), dst_index(0),
            periodic_shift(IntVect(AMREX_D_DECL(0, 0, 0)))
            {}

        bool operator<(const NeighborCopyTag& other) const {
            if (level != other.level) return level < other.level;
            if (grid != other.grid) return grid < other.grid;
            if (tile != other.tile) return tile < other.tile;
            AMREX_D_TERM(
            if (periodic_shift[0] != other.periodic_shift[0])
                return periodic_shift[0] < other.periodic_shift[0];,
            if (periodic_shift[1] != other.periodic_shift[1])
                return periodic_shift[1] < other.periodic_shift[1];,
            if (periodic_shift[2] != other.periodic_shift[2])
                return periodic_shift[2] < other.periodic_shift[2];
                )
            return false;
        }

        bool operator==(const NeighborCopyTag& other) const {
            return (level == other.level) and (grid == other.grid) and (tile == other.tile)
                AMREX_D_TERM(
                    and (periodic_shift[0] == other.periodic_shift[0]),
                    and (periodic_shift[1] == other.periodic_shift[1]),
                    and (periodic_shift[2] == other.periodic_shift[2])
                    );
        }

        bool operator!=(const NeighborCopyTag& other) const
        {
            return !operator==(other);
        }

        friend std::ostream& operator<< (std::ostream& os, const NeighborCopyTag& tag)
        {
            os << tag.level << " " << tag.grid << " " << tag.tile << " " << tag.periodic_shift;
            if (!os.good())
                amrex::Error("operator<<(ostream&,Particle<NReal, NInt>&) failed");
            return os;
        }
    };

    struct NeighborCommTag {

        NeighborCommTag(int pid, int lid, int gid, int tid)
            : proc_id(pid), level_id(lid), grid_id(gid), tile_id(tid)
            {}

        int proc_id;
        int level_id;
        int grid_id;
        int tile_id;

        bool operator<(const NeighborCommTag& other) const {
            return (proc_id < other.proc_id ||
                    (proc_id == other.proc_id &&
                     grid_id < other.grid_id) ||
                    (proc_id == other.proc_id &&
                     grid_id == other.grid_id &&
                     tile_id < other.tile_id) ||
                    (proc_id == other.proc_id &&
                     grid_id == other.grid_id &&
                     tile_id == other.tile_id &&
                     level_id < other.level_id));
        }

        bool operator==(const NeighborCommTag& other) const {
            return ( (proc_id == other.proc_id) and
                     (grid_id == other.grid_id) and
                     (tile_id == other.tile_id) and
                     (level_id == other.level_id));
        }

        friend std::ostream& operator<< (std::ostream& os, const NeighborCommTag& tag)
        {
            os << tag.proc_id << " " << tag.level_id << " " << tag.grid_id << " " << tag.tile_id;
            if (!os.good())
                amrex::Error("operator<<(ostream&,Particle<NReal, NInt>&) failed");
            return os;
        }
    };

public:

    using ParticleType = typename ParticleContainer<NStructReal,
                                                    NStructInt, 0, 0>::ParticleType;
    using MyParIter = ParIter<NStructReal, NStructInt, 0, 0>;
    using PairIndex = std::pair<int, int>;
    using NeighborCommMap = std::map<NeighborCommTag, Vector<char> >;
    using AoS = typename ParticleContainer<NStructReal, NStructInt, 0, 0>::AoS;
    using ParticleVector = typename ParticleContainer<NStructReal, NStructInt, 0, 0>::ParticleVector;
    using IntVector  = typename ParticleContainer<NStructReal, NStructInt, 0, 0>::IntVector;
    using SendBuffer = typename ParticleContainer<NStructReal, NStructInt, 0, 0>::SendBuffer;

    NeighborParticleContainer(ParGDBBase* gdb, int ncells);

    NeighborParticleContainer(const Geometry            & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba,
                              int                         nneighbor);

    NeighborParticleContainer(const Vector<Geometry>            & geom,
                              const Vector<DistributionMapping> & dmap,
                              const Vector<BoxArray>            & ba,
                              const Vector<int>                 & rr,
                              int                               nneighbor);

    ///
    /// Regrid functions
    ///
    void Regrid(const DistributionMapping& dmap, const BoxArray& ba);
    void Regrid(const DistributionMapping& dmap, const BoxArray& ba, const int lev);
    void Regrid(const Vector<DistributionMapping>& dmap, const Vector<BoxArray>& ba);

    ///
    /// This fills the neighbor buffers for each tile with the proper data
    ///
    void fillNeighbors();

    ///
    /// This updates the neighbors with their current particle data.
    ///
    void updateNeighbors();

    ///
    /// Each tile clears its neighbors, freeing the memory
    ///
    void clearNeighbors();

    ///
    /// Build a Neighbor List for each tile
    ///
    template <class CheckPair>
    void buildNeighborList(CheckPair check_pair, bool sort=false);

    void printNeighborList();

    void setRealCommComp(int i, bool value);
    void setIntCommComp (int i, bool value);

    ParticleVector& GetNeighbors(int lev, int grid, int tile)
    {
        return neighbors[lev][std::make_pair(grid,tile)];
    }

    const ParticleVector& GetNeighbors(int lev, int grid, int tile) const
    {
        return neighbors[lev][std::make_pair(grid,tile)];
    }

    IntVector& GetNeighborList(int lev, int grid, int tile)
    {
        return neighbor_list[lev][std::make_pair(grid,tile)];
    }

    const IntVector& GetNeighborList(int lev, int grid, int tile) const
    {
        return neighbor_list[lev][std::make_pair(grid,tile)];
    }

    void RedistributeLocal()
    {
        const int lev_min = 0;
        const int lev_max = 0;
        const int nGrow = 0;
        const int local = 1;
        clearNeighbors();
        this->Redistribute(lev_min, lev_max, nGrow, local);    
    }

#ifdef AMREX_USE_CUDA

    template <class CheckPair>
    void buildNeighborListGPU (CheckPair check_pair);
    
    void printNeighborListGPU ();
    
    void sortParticlesByNeighborDestGPU ();
    
    void fillNeighborsGPU ();
    
    void updateNeighborsGPU ();
    
    void clearNeighborsGPU ();
    
#else 

    template <class CheckPair>
    void buildNeighborListCPU (CheckPair check_pair, bool sort);
    
    void printNeighborListCPU(const std::string& prefix);
    
    void fillNeighborsCPU ();
    
    void updateNeighborsCPU (bool reuse_rcv_counts=true);
    
    void clearNeighborsCPU ();

#endif

protected:

    void cacheNeighborInfo();

    ///
    /// This builds the internal mask data structure used for looking up neighbors
    ///
    void BuildMasks();

    void GetNeighborCommTags();

    void GetCommTagsBox(Vector<NeighborCommTag>& tags, const int lev, const Box& in_box);

    void resizeContainers(const int lev);

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

    void GetNeighborCommTags(Vector<NeighborCommTag>& tags, const int lev, Box box);

    void getNeighborTags(Vector<NeighborCopyTag>& tags,
                         const ParticleType& p,
                         const int nGrow, const NeighborCopyTag& src_tag, const MyParIter& pti);

    void getNeighborTags(Vector<NeighborCopyTag>& tags,
                         const ParticleType& p,
                         const IntVect& nGrow, const NeighborCopyTag& src_tag, const MyParIter& pti);

    IntVect computeRefFac(const int src_lev, const int lev);

    amrex::Vector<std::map<PairIndex, ParticleVector> > neighbors;
    amrex::Vector<std::map<PairIndex, IntVector> >      neighbor_list;
    const size_t pdata_size = sizeof(ParticleType);

    static constexpr int num_mask_comps = 3;  //!< grid, tile, level
    size_t cdata_size;
    int m_num_neighbor_cells;
    amrex::Vector<NeighborCommTag> local_neighbors;
    amrex::Vector<std::unique_ptr<iMultiFab> > mask_ptr;

    amrex::Vector<std::map<PairIndex, Vector<Vector<NeighborCopyTag> > > > buffer_tag_cache;
    amrex::Vector<std::map<PairIndex, int> > local_neighbor_sizes;

    //! each proc knows how many sends it will do, and how many bytes it will rcv
    //! from each other proc.
    amrex::Vector<int> neighbor_procs;
    amrex::Vector<long> rcvs;
    long num_snds;
    std::map<int, Vector<char> > send_data;

    std::array<bool, AMREX_SPACEDIM + NStructReal> rc;
    std::array<bool, 2 + NStructInt>  ic;

    static bool use_mask;

#ifdef AMREX_USE_CUDA
    void fillNeighborsMPIGPU(std::map<int, SendBuffer>& not_ours);
    
    void buildNeighborMask ();

    struct NeighborCode {
        int grid_id;
        IntVect periodic_shift;
    };

    struct NeighborTask {
        int grid_id;
        Box box;
        IntVect periodic_shift;

        NeighborTask(int a_grid_id, const Box& a_box, const IntVect& a_periodic_shift)
            : grid_id(a_grid_id), box(a_box), periodic_shift(a_periodic_shift) {}

        bool operator<(const NeighborTask& other) const {
            if (grid_id != other.grid_id) return grid_id < other.grid_id;
            if (box     != other.box    ) return box < other.box;
            AMREX_D_TERM(
            if (periodic_shift[0] != other.periodic_shift[0])
                return periodic_shift[0] < other.periodic_shift[0];,
            if (periodic_shift[1] != other.periodic_shift[1])
                return periodic_shift[1] < other.periodic_shift[1];,
            if (periodic_shift[2] != other.periodic_shift[2])
                return periodic_shift[2] < other.periodic_shift[2];
                )
            return false;
        }
    };
        
    // These are used to keep track of which particles need to be ghosted to which grids
    bool m_neighbor_mask_initialized = false;
    std::unique_ptr<amrex::iMultiFab> m_neighbor_mask_ptr;
    std::map<int, std::vector<std::vector<NeighborCode> > > m_grid_map;
    std::map<int, amrex::Gpu::HostVector<int> > m_start;
    std::map<int, amrex::Gpu::HostVector<int> > m_stop;
    
    std::map<std::pair<int, int>, amrex::NeighborList<ParticleType> > m_neighbor_list;
#endif
};
    
#include "AMReX_NeighborParticlesI.H"
    
#ifdef AMREX_USE_CUDA
#include "AMReX_NeighborParticlesGPUImpl.H"
#else
#include "AMReX_NeighborParticlesCPUImpl.H"
#endif
    
}

#endif // _NEIGHBORPARTICLES_H_
