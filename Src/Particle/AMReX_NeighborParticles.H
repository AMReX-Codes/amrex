#ifndef AMREX_NEIGHBORPARTICLES_H_
#define AMREX_NEIGHBORPARTICLES_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_ParticleUtil.H>
#include <AMReX_NeighborList.H>
#include <AMReX_OpenMP.H>
#include <AMReX_ParticleTile.H>

namespace amrex {

  struct NeighborCode
  {
      int grid_id;
      IntVect periodic_shift;
  };

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
template <int NStructReal, int NStructInt, int NArrayReal=0, int NArrayInt=0>
class NeighborParticleContainer // NOLINT(cppcoreguidelines-virtual-class-destructor) // clang-tidy seems wrong.
    : public ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
{
public:
    using ParticleContainerType = ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    using ParticleType = typename ParticleContainerType::ParticleType;
    using SuperParticleType = typename ParticleContainerType::SuperParticleType;
    using NeighborListContainerType = Vector<std::map<std::pair<int, int>, amrex::NeighborList<ParticleType> > >;
private:
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
                amrex::Error("operator<<(ostream&, const NeighborIndexMap& nim) failed");
            return os;
        }

    };

    struct NeighborCopyTag {
        int level = -1;
        int grid = -1;
        int tile = -1;
        int src_index = 0;
        int dst_index = 0;
        IntVect periodic_shift = IntVect(0);

        NeighborCopyTag () = default;

        NeighborCopyTag (int a_level, int a_grid, int a_tile) :
            level(a_level), grid(a_grid), tile(a_tile)
            {}

        bool operator< (const NeighborCopyTag& other) const {
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

        bool operator== (const NeighborCopyTag& other) const {
            return (level == other.level) && (grid == other.grid) && (tile == other.tile)
                AMREX_D_TERM(
                    && (periodic_shift[0] == other.periodic_shift[0]),
                    && (periodic_shift[1] == other.periodic_shift[1]),
                    && (periodic_shift[2] == other.periodic_shift[2])
                    );
        }

        bool operator!= (const NeighborCopyTag& other) const
        {
            return !operator==(other);
        }

        friend std::ostream& operator<< (std::ostream& os, const NeighborCopyTag& tag)
        {
            os << tag.level << " " << tag.grid << " " << tag.tile << " " << tag.periodic_shift;
            if (!os.good())
                amrex::Error("operator<<(ostream&, const NeighborCopyTag&) failed");
            return os;
        }
    };

    struct InverseCopyTag
    {
        int src_grid;
        int src_tile;
        int src_index;
        int src_level;
    };

    friend std::ostream& operator<< (std::ostream& os, const InverseCopyTag& tag)
    {
        os << tag.src_level << " " << tag.src_grid << " " << tag.src_tile << " " << tag.src_index;
        if (!os.good())
            amrex::Error("operator<<(ostream&, const InverseCopyTag&) failed");
        return os;
    }

    struct NeighborCommTag {

        NeighborCommTag (int pid, int lid, int gid, int tid)
            : proc_id(pid), level_id(lid), grid_id(gid), tile_id(tid)
            {}

        int proc_id;
        int level_id;
        int grid_id;
        int tile_id;

        bool operator< (const NeighborCommTag& other) const {
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

        bool operator== (const NeighborCommTag& other) const {
            return ( (proc_id == other.proc_id) &&
                     (grid_id == other.grid_id) &&
                     (tile_id == other.tile_id) &&
                     (level_id == other.level_id));
        }

        friend std::ostream& operator<< (std::ostream& os, const NeighborCommTag& tag)
        {
            os << tag.proc_id << " " << tag.level_id << " " << tag.grid_id << " " << tag.tile_id;
            if (!os.good())
                amrex::Error("operator<<(ostream&, const NeighborCommTag&) failed");
            return os;
        }
    };

public:

    using MyParIter = ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    using PairIndex = std::pair<int, int>;
    using NeighborCommMap = std::map<NeighborCommTag, Vector<char> >;
    using AoS = typename ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::AoS;
    using ParticleVector = typename ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::ParticleVector;
    using ParticleTile = typename ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::ParticleTileType;
    using IntVector  = typename ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::IntVector;
    using SendBuffer = typename ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::SendBuffer;

    NeighborParticleContainer (ParGDBBase* gdb, int ncells);

    NeighborParticleContainer (const Geometry            & geom,
                               const DistributionMapping & dmap,
                               const BoxArray            & ba,
                               int                         nneighbor);

    NeighborParticleContainer (const Vector<Geometry>            & geom,
                               const Vector<DistributionMapping> & dmap,
                               const Vector<BoxArray>            & ba,
                               const Vector<int>                 & rr,
                               int                               nneighbor);

    ~NeighborParticleContainer () override = default;

    NeighborParticleContainer ( const NeighborParticleContainer &) = delete;
    NeighborParticleContainer& operator= ( const NeighborParticleContainer & ) = delete;

    NeighborParticleContainer ( NeighborParticleContainer && ) = default; // NOLINT(performance-noexcept-move-constructor)
    NeighborParticleContainer& operator= ( NeighborParticleContainer && ) = default; // NOLINT(performance-noexcept-move-constructor)

    ///
    /// Regrid functions
    ///
    void Regrid (const DistributionMapping& dmap, const BoxArray& ba);
    void Regrid (const DistributionMapping& dmap, const BoxArray& ba, int lev);
    void Regrid (const Vector<DistributionMapping>& dmap, const Vector<BoxArray>& ba);

    ///
    /// This fills the neighbor buffers for each tile with the proper data
    ///
    void fillNeighbors ();

    ///
    /// This does an "inverse" fillNeighbors operation, meaning that it adds
    /// data from the ghost particles to the corresponding real ones.
    ///
    void sumNeighbors (int real_start_comp, int real_num_comp,
                       int int_start_comp, int int_num_comp);

    ///
    /// This updates the neighbors with their current particle data.
    ///
    void updateNeighbors (bool boundary_neighbors_only=false);

    ///
    /// Each tile clears its neighbors, freeing the memory
    ///
    void clearNeighbors ();

    ///
    /// Build a Neighbor List for each tile
    ///
    template <class CheckPair>
    void buildNeighborList (CheckPair&& check_pair, bool sort=false);

    ///
    /// Build a Neighbor List for each tile
    ///
    template <class CheckPair, class OtherPCType>
    void buildNeighborList (CheckPair&& check_pair, OtherPCType& other,
                            Vector<std::map<std::pair<int, int>, amrex::NeighborList<typename OtherPCType::ParticleType> > >& neighbor_lists,
                            bool sort=false);

    ///
    /// Build a Neighbor List for each tile
    ///
    template <class CheckPair>
    void buildNeighborList (CheckPair&& check_pair, int type_ind, int* ref_ratio,
                            int num_bin_types=1, bool sort=false);

    template <class CheckPair>
    void selectActualNeighbors (CheckPair&& check_pair, int num_cells=1);

    void printNeighborList ();

    void setRealCommComp (int i, bool value);
    void setIntCommComp (int i, bool value);

    ParticleTile& GetNeighbors (int lev, int grid, int tile)
    {
        return neighbors[lev][std::make_pair(grid,tile)];
    }

    const ParticleTile& GetNeighbors (int lev, int grid, int tile) const
    {
        return neighbors[lev][std::make_pair(grid,tile)];
    }

    template <typename T,
              typename std::enable_if<std::is_same<T,bool>::value,int>::type=0>
    void AddRealComp (T communicate=true)
    {
        ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
            AddRealComp(communicate);
        ghost_real_comp.push_back(communicate);
        calcCommSize();
    }

    template <typename T,
              typename std::enable_if<std::is_same<T,bool>::value,int>::type=0>
    void AddIntComp (T communicate=true)
    {
        ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
            AddIntComp(communicate);
        ghost_int_comp.push_back(communicate);
        calcCommSize();
    }

    void Redistribute (int lev_min=0, int lev_max=-1, int nGrow=0, int local=0)
    {
        clearNeighbors();
        ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
            ::Redistribute(lev_min, lev_max, nGrow, local);
    }

    void RedistributeLocal ()
    {
        const int lev_min = 0;
        const int lev_max = 0;
        const int nGrow = 0;
        const int local = 1;
        clearNeighbors();
        this->Redistribute(lev_min, lev_max, nGrow, local);
    }

#ifdef AMREX_USE_GPU
    void fillNeighborsGPU ();
    void updateNeighborsGPU (bool boundary_neighbors_only=false);
    void clearNeighborsGPU ();
#else
    void fillNeighborsCPU ();
    void sumNeighborsCPU (int real_start_comp, int real_num_comp,
                          int int_start_comp, int int_num_comp);
    void updateNeighborsCPU (bool reuse_rcv_counts=true);
    void clearNeighborsCPU ();
#endif

    void setEnableInverse (bool flag)
    {
        enable_inverse = flag;
        calcCommSize();
    }

    bool enableInverse () { return enable_inverse; }

    void buildNeighborMask ();

    void buildNeighborCopyOp (bool use_boundary_neighbor=false);

protected:

    void cacheNeighborInfo ();

    ///
    /// This builds the internal mask data structure used for looking up neighbors
    ///
    void BuildMasks ();

    ///
    /// Are the masks computed by the above function still valid?
    ///
    bool areMasksValid ();

    void GetNeighborCommTags ();

    void GetCommTagsBox (Vector<NeighborCommTag>& tags, int lev, const Box& in_box);

    void resizeContainers (int num_levels);

    void initializeCommComps ();

    void calcCommSize ();

    ///
    /// Perform the MPI communication necessary to fill neighbor buffers
    ///
    void fillNeighborsMPI (bool reuse_rcv_counts);

    void sumNeighborsMPI (std::map<int, Vector<char> >& not_ours,
                          int real_start_comp, int real_num_comp,
                          int int_start_comp, int int_num_comp);

    ///
    /// Perform handshake to figure out how many bytes each proc should receive
    ///
    void getRcvCountsMPI ();

    void getNeighborTags (Vector<NeighborCopyTag>& tags,
                         const ParticleType& p,
                          int nGrow, const NeighborCopyTag& src_tag, const MyParIter& pti);

    void getNeighborTags (Vector<NeighborCopyTag>& tags,
                         const ParticleType& p,
                         const IntVect& nGrow, const NeighborCopyTag& src_tag, const MyParIter& pti);

    IntVect computeRefFac (int src_lev, int lev);

    Vector<std::map<PairIndex, Vector<InverseCopyTag> > > inverse_tags;
    Vector<std::map<PairIndex, ParticleTile> > neighbors;
    Vector<std::map<PairIndex, IntVector> >      neighbor_list;
    const size_t pdata_size = sizeof(ParticleType);

    static constexpr int num_mask_comps = 3;  //!< grid, tile, level
    size_t cdata_size;
    int m_num_neighbor_cells;
    Vector<NeighborCommTag> local_neighbors;
    Vector<std::unique_ptr<iMultiFab> > mask_ptr;

    Vector<std::map<PairIndex, Vector<Vector<NeighborCopyTag> > > > buffer_tag_cache;
    Vector<std::map<PairIndex, int> > local_neighbor_sizes;

    //! each proc knows how many sends it will do, and how many bytes it will rcv
    //! from each other proc.
    Vector<int> neighbor_procs;
    Vector<Long> rcvs;
    Long num_snds;
    std::map<int, Vector<char> > send_data;

    Vector<int> ghost_real_comp;
    Vector<int> ghost_int_comp;

    static bool use_mask;

    static bool enable_inverse;

#ifdef AMREX_USE_GPU

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
    // {pboxid: [[neighbor codes] size= #intersected nbrs] size= #boundary boxes}
    std::map<int, std::vector<std::vector<NeighborCode> > > m_grid_map;

    // {pboxid: [ith intersected pbox: nbr codes]}
    std::map<int, amrex::Gpu::DeviceVector<NeighborCode> > m_code_array;
    // {pboxid: [ith intersected pbox: intersection]}
    std::map<int, amrex::Gpu::DeviceVector<Box> >          m_isec_boxes;
    std::map<int, amrex::Gpu::DeviceVector<int> > m_code_offsets; // to be removed

    ParticleCopyOp neighbor_copy_op;
    ParticleCopyPlan neighbor_copy_plan;

    amrex::PODVector<char, PolymorphicArenaAllocator<char> > snd_buffer;
    amrex::Gpu::DeviceVector<char> rcv_buffer;

    Gpu::PinnedVector<char> pinned_snd_buffer;
    Gpu::PinnedVector<char> pinned_rcv_buffer;
#endif

    NeighborListContainerType m_neighbor_list;

    Vector<std::map<std::pair<int, int>, amrex::Gpu::DeviceVector<int> > > m_boundary_particle_ids;

    [[nodiscard]] bool hasNeighbors() const { return m_has_neighbors; }

    bool m_has_neighbors = false;
};

#include "AMReX_NeighborParticlesI.H"

#ifdef AMREX_USE_GPU
#include "AMReX_NeighborParticlesGPUImpl.H"
#else
#include "AMReX_NeighborParticlesCPUImpl.H"
#endif

}

#endif // _NEIGHBORPARTICLES_H_
