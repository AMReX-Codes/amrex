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
#include <AMReX_ParticleTileRT.H>

namespace amrex {

  struct NeighborCode
  {
      int grid_id;
      Box grid_box;
      IntVect periodic_shift;
  };

namespace detail {

/// Component lists for the same-rank ghost refresh on the CPU
/// (updateNeighborsCPU). Only the components flagged for ghost communication
/// are copied, and from a group of per-body "slot" arrays (contact history)
/// only the occupied slots. A slot is empty when both its keys are -1; the key
/// arrays and the payload arrays are runtime components, the payload being
/// nblocks blocks of nslots entries each (block b, slot s at real0 + b*nslots + s).
struct GhostCopyPlan
{
    Vector<int> ct_real;   ///< flagged compile-time array reals (index into rdata)
    Vector<int> rt_real;   ///< flagged runtime reals outside the slot group
    Vector<int> ct_int;    ///< flagged compile-time array ints (index into idata)
    Vector<int> rt_int;    ///< flagged runtime ints outside the slot group
    int nslots  = 0;       ///< slots per body, 0 = no slot group
    int key_cpu = 0;       ///< first runtime int of the slot cpu keys
    int key_id  = 0;       ///< first runtime int of the slot id keys
    int real0   = 0;       ///< first runtime real of the slot payload
    int nblocks = 0;       ///< payload blocks of nslots reals

    /// Layout of a ghost record on the MPI path, filled only when the masked MPI
    /// path is on (setGhostMpiMasked) and a slot group is declared: packing and
    /// unpacking then skip the payload of the empty slots (the record keeps its
    /// size and layout; the bytes of an empty slot's payload are neither written
    /// by the sender nor read by the receiver).
    bool wire = false;
    Vector<int> wire_real_slot;     ///< per real component: slot of a payload component, else -1
    Vector<int> wire_key_cpu_off;   ///< per slot: byte offset of the cpu key inside a record
    Vector<int> wire_key_id_off;    ///< per slot: byte offset of the id key inside a record

    /// Batched packing and unpacking of the records (setGhostBatched): the byte
    /// offset of every component inside a record, so that one component at a time
    /// can be moved for all the ghosts of a message. The records themselves are
    /// the ones of the masked MPI path. Filled together with the layout above for
    /// particles stored as an array of structures without extra compile-time
    /// components, when the record size is the communication size.
    bool batch = false;
    int  rec_size = 0;                       ///< bytes of a record
    int  wire_pos_off[3] = {-1, -1, -1};     ///< array-of-structures position components
    int  wire_idcpu_off[2] = {-1, -1};       ///< id and cpu parts of the packed id/cpu word
    Vector<int> wire_real_off;               ///< per real component: byte offset in a record, -1 if not sent
    Vector<int> wire_int_off;                ///< per int component: byte offset in a record, -1 if not sent

    /// Compact records (setGhostCompact), available when the batched layout above is:
    /// the head of a ghost holds every component that is not the payload of a slot, in
    /// the order of the record, and at its end the offset (an int) of the payload of the
    /// ghost in the payload area at the end of the message. The payload of the occupied
    /// slots lies there, cpay_bytes (nblocks reals) each, in slot order. A slot is
    /// occupied when its two keys, which are in the head, are not both -1.
    bool compact = false;
    int  chead_size = 0;                     ///< bytes of a head
    int  cpoff_off  = -1;                    ///< byte offset of the payload offset in a head
    int  cpay_bytes = 0;                     ///< bytes of the payload of one occupied slot
    int  cwire_pos_off[3] = {-1, -1, -1};    ///< array-of-structures position components
    int  cwire_idcpu_off[2] = {-1, -1};      ///< id and cpu parts of the packed id/cpu word
    Vector<int> cwire_real_off;              ///< per real component: byte offset in a head, -1 if not sent or in the payload
    Vector<int> cwire_int_off;               ///< per int component: byte offset in a head, -1 if not sent
    Vector<int> cwire_key_cpu_off;           ///< per slot: byte offset of the cpu key in a head
    Vector<int> cwire_key_id_off;            ///< per slot: byte offset of the id key in a head
};

/// Largest slot group handled by the masked MPI path (the occupancy of a record
/// is kept in a stack array of this size).
inline constexpr int MaxWireSlots = 1024;

/// Most payload blocks (reals per slot) handled by the compact records.
inline constexpr int MaxPayloadBlocks = 32;

/// Payload area of the message to one rank (compact records): where it starts
/// and the bytes of it used by the ghosts packed so far in this update.
struct GhostSendInfo
{
    char* pay = nullptr;
    std::size_t used = 0;
};

/// Work arrays of the batched packing and unpacking, kept between calls.
struct GhostBatchScratch
{
    Vector<char*>   rec;     ///< pack: record of each remote ghost of the source tile
    Vector<int>     si;      ///< pack: its source index
    Vector<IntVect> shift;   ///< pack: its periodic shift
    Vector<int>     lists;   ///< per slot, the ghosts of the batch where the slot is occupied (nslots lists of n entries)
    Vector<int>     cnt;     ///< length of each list
    // compact records
    Vector<GhostSendInfo*> sinfo; ///< pack: payload area of the message of each remote ghost
    Vector<int>     ord;     ///< occupied slots counted so far, per ghost of the batch
    Vector<int>     kord;    ///< like lists: the place of the slot among the occupied slots of the ghost
    Vector<int>     poff;    ///< payload offset of each ghost of the batch
};

} // namespace detail

///
/// This is a container for particles that undergo short-range interactions.
/// In addition to the normal ParticleContainer methods, each tile contains a "neighbor
/// buffer" that is filled with data corresponding to the particles within 1 neighbor cell of
/// the tile boundaries. This allows the N^2 search over each pair of particles to proceed
/// locally, instead of over the entire domain.
///
/// Note that neighbor particles are different than "ghost" particles, which are used
/// in AMR subcycling to keep track of coarse level particles that may move on to fine
/// levels during a fine level time step.
///
template <typename T_ParticleType, int T_NArrayReal=0, int T_NArrayInt=0,
          template<class> class Allocator=DefaultAllocator,
          class CellAssignor=DefaultAssignor>
class NeighborParticleContainer_impl // NOLINT(cppcoreguidelines-virtual-class-destructor) // clang-tidy seems wrong.
    : public ParticleContainer_impl<T_ParticleType, T_NArrayReal, T_NArrayInt, Allocator, CellAssignor>
{
public:
    using ParticleContainerType =
        ParticleContainer_impl<T_ParticleType, T_NArrayReal, T_NArrayInt, Allocator, CellAssignor>;
    using ParticleType = typename ParticleContainerType::ParticleType;
    using SuperParticleType = typename ParticleContainerType::SuperParticleType;
    using NeighborListContainerType =
        Vector<std::map<std::pair<int, int>,
            amrex::NeighborList<ParticleType> > >;
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
            if (!os.good()) {
                amrex::Error("operator<<(ostream&, const NeighborIndexMap& nim) failed");
            }
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
            if (level != other.level) { return level < other.level; }
            if (grid != other.grid) { return grid < other.grid; }
            if (tile != other.tile) { return tile < other.tile; }
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
            if (!os.good()) {
                amrex::Error("operator<<(ostream&, const NeighborCopyTag&) failed");
            }
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
        if (!os.good()) {
            amrex::Error("operator<<(ostream&, const InverseCopyTag&) failed");
        }
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
            if (!os.good()) {
                amrex::Error("operator<<(ostream&, const NeighborCommTag&) failed");
            }
            return os;
        }
    };

public:

    static constexpr int NArrayReal = ParticleContainerType::NArrayReal;
    static constexpr int NArrayInt = ParticleContainerType::NArrayInt;
    static constexpr int NStructReal = ParticleContainerType::NStructReal;
    static constexpr int NStructInt = ParticleContainerType::NStructInt;
    static constexpr int RealCommCompStart = ParticleType::is_soa_particle ? 0 : AMREX_SPACEDIM + NStructReal;
    static constexpr int IntCommCompStart = 2 + NStructInt;

    using MyParIter = typename ParticleContainerType::ParIterType;
    using PairIndex = std::pair<int, int>;
    using NeighborCommMap = std::map<NeighborCommTag, Vector<char> >;
    using AoS = typename ParticleContainerType::AoS;
    using ParticleVector = typename ParticleContainerType::ParticleVector;
    using ParticleTile = typename ParticleContainerType::ParticleTileType;
    using IntVector  = typename ParticleContainerType::IntVector;
    using PTDType = typename ParticleContainerType::PTDType;
    using ConstPTDType = typename ParticleContainerType::ConstPTDType;

    NeighborParticleContainer_impl (ParGDBBase* gdb, int ncells);

    NeighborParticleContainer_impl (const Geometry            & geom,
                                    const DistributionMapping & dmap,
                                    const BoxArray            & ba,
                                    int                         nneighbor);

    NeighborParticleContainer_impl (const Vector<Geometry>            & geom,
                                    const Vector<DistributionMapping> & dmap,
                                    const Vector<BoxArray>            & ba,
                                    const Vector<int>                 & rr,
                                    int                               nneighbor);

    ~NeighborParticleContainer_impl () override = default;

    NeighborParticleContainer_impl ( const NeighborParticleContainer_impl &) = delete;
    NeighborParticleContainer_impl& operator= ( const NeighborParticleContainer_impl & ) = delete;

    NeighborParticleContainer_impl ( NeighborParticleContainer_impl && ) = default; // NOLINT(performance-noexcept-move-constructor)
    NeighborParticleContainer_impl& operator= ( NeighborParticleContainer_impl && ) = default; // NOLINT(performance-noexcept-move-constructor)

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
    /// Fill neighbor buffers, but only communicate particles within
    /// physical distance search_radius of tile boundaries.
    /// Reduces communication volume when search_radius << cell size.
    /// CPU build only; the GPU path is not yet implemented.
    ///
    void fillNeighbors (amrex::Real search_radius);

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
    void buildNeighborList (CheckPair const& check_pair, bool sort=false);

    ///
    /// Build a Neighbor List for each tile, using custom bin_size
    /// instead of the mesh cell size. Useful when the interaction
    /// radius is much smaller than dx.
    ///
    template <class CheckPair>
    void buildNeighborList (CheckPair const& check_pair, amrex::Real bin_size, bool sort=false);

    ///
    /// Build a Neighbor List for each tile
    ///
    template <class CheckPair, class OtherPCType>
    void buildNeighborList (CheckPair const& check_pair, OtherPCType& other,
                            Vector<std::map<std::pair<int, int>,
                                amrex::NeighborList<typename OtherPCType::ParticleType> > >& neighbor_lists,
                            bool sort=false);

    ///
    /// Build a cross-container Neighbor List for each tile, using custom
    /// bin_size instead of the mesh cell size. Useful when the interaction
    /// radius is much smaller than dx.
    ///
    template <class CheckPair, class OtherPCType>
    void buildNeighborList (CheckPair const& check_pair, OtherPCType& other,
                            Vector<std::map<std::pair<int, int>,
                                amrex::NeighborList<typename OtherPCType::ParticleType> > >& neighbor_lists,
                            amrex::Real bin_size, bool sort=false);

    ///
    /// Build a Neighbor List for each tile
    ///
    template <class CheckPair>
    void buildNeighborList (CheckPair const& check_pair, int type_ind, int* ref_ratio,
                            int num_bin_types=1, bool sort=false);

    template <class CheckPair>
    void selectActualNeighbors (CheckPair const& check_pair, int num_cells=1);

    void printNeighborList ();

    void setRealCommComp (int i, bool value);
    void setIntCommComp (int i, bool value);

    ///
    /// Declare a group of per-body slot arrays among the runtime components
    /// (for example contact history): nslots slots per body, the cpu keys at
    /// runtime int key_cpu + s, the id keys at key_id + s, and nblocks blocks
    /// of nslots reals from runtime real real0 (block b, slot s at
    /// real0 + b*nslots + s). A slot is empty when both its keys are -1. With
    /// the group declared, the masked same-rank ghost refresh on the CPU (see
    /// setGhostCopyMasked) copies the payload of the occupied slots only. All
    /// the components of the group must be flagged for ghost communication,
    /// otherwise the group is ignored.
    ///
    void setGhostSlotGroup (int nslots, int key_cpu, int key_id, int real0, int nblocks)
    {
        ghost_slot_n       = nslots;
        ghost_slot_key_cpu = key_cpu;
        ghost_slot_key_id  = key_id;
        ghost_slot_real0   = real0;
        ghost_slot_nblocks = nblocks;
    }

    ///
    /// true: the same-rank ghost refresh on the CPU (updateNeighbors) copies
    /// only the components flagged for ghost communication, as the MPI path
    /// does, and from a declared slot group only the occupied slots.
    /// false (default): it copies every component of every ghost.
    ///
    void setGhostCopyMasked (bool flag) { ghost_copy_masked = flag; }
    [[nodiscard]] bool ghostCopyMasked () const { return ghost_copy_masked; }

    ///
    /// true: with the masked copy on and a slot group declared, the packing and
    /// unpacking of the ghost records on the MPI path (CPU) skip the payload of
    /// the empty slots. The records keep their size and layout.
    /// false (default): every component of a record is packed and unpacked.
    ///
    void setGhostMpiMasked (bool flag) { ghost_mpi_masked = flag; }
    [[nodiscard]] bool ghostMpiMasked () const { return ghost_mpi_masked; }

    ///
    /// true (CPU): the receive buffer of the ghost exchange is kept between
    /// updates instead of being allocated and cleared for each one, and the
    /// component lists of the masked copy are made once per rebuild instead of
    /// once per update. The ghosts do not change.
    /// false (default): as before.
    ///
    void setGhostPersistentRecv (bool flag) { ghost_persistent_recv = flag; ghost_plan_valid = false; }
    [[nodiscard]] bool ghostPersistentRecv () const { return ghost_persistent_recv; }

    ///
    /// true (CPU): an update of the ghosts without rebuild (updateNeighbors)
    /// writes them straight into the particle tiles, without the staging tiles
    /// and the final copy. Not with enableInverse(). The ghosts do not change.
    /// false (default): as before.
    ///
    void setGhostDirect (bool flag) { ghost_direct = flag; }
    [[nodiscard]] bool ghostDirect () const { return ghost_direct; }

    ///
    /// true (CPU): with the masked MPI path on and a slot group declared, the
    /// records of a message are packed and unpacked one component at a time for
    /// all of its ghosts, for particles stored as an array of structures without
    /// extra compile-time components. The records do not change.
    /// false (default): one ghost at a time.
    ///
    void setGhostBatched (bool flag) { ghost_batched = flag; ghost_plan_valid = false; }
    [[nodiscard]] bool ghostBatched () const { return ghost_batched; }

    ///
    /// true (CPU): the sends of the ghost exchange are posted without blocking
    /// and waited for together with the receives. The data do not change.
    /// false (default): one blocking send after another.
    ///
    void setGhostIsend (bool flag) { ghost_isend = flag; }
    [[nodiscard]] bool ghostIsend () const { return ghost_isend; }

    ///
    /// true (CPU): with the masked MPI path and the batched packing on, a slot
    /// group declared, one thread per rank and no inverse tags, a ghost travels as a
    /// fixed-size head (every component but the payload of the slots, and the
    /// place of its payload) and the payload of its occupied slots, which is
    /// appended to an area at the end of the message; only the used part of a
    /// message is sent. The ghosts do not change. Takes effect at the next rebuild
    /// of the neighbor lists.
    /// false (default): fixed-size records.
    ///
    void setGhostCompact (bool flag) { ghost_compact = flag; ghost_plan_valid = false; }
    [[nodiscard]] bool ghostCompact () const { return ghost_compact; }

    ParticleTile& GetNeighbors (int lev, int grid, int tile)
    {
        return neighbors[lev][std::make_pair(grid,tile)];
    }

    const ParticleTile& GetNeighbors (int lev, int grid, int tile) const
    {
        return neighbors[lev][std::make_pair(grid,tile)];
    }

    template <std::same_as<bool> T>
    void AddRealComp (T communicate=true)
    {
        ParticleContainerType::AddRealComp(communicate);
        ghost_real_comp.push_back(communicate);
        calcCommSize();
    }

    template <std::same_as<bool> T>
    void AddIntComp (T communicate=true)
    {
        ParticleContainerType::AddIntComp(communicate);
        ghost_int_comp.push_back(communicate);
        calcCommSize();
    }

    void Redistribute (int lev_min=0, int lev_max=-1, int nGrow=0,
                       int local=0, bool remove_negative=true)
    {
        clearNeighbors();
        ParticleContainerType::Redistribute(lev_min, lev_max, nGrow, local, remove_negative);
    }

    void Redistribute (int lev_min, int lev_max, IntVect nGrow,
                       bool local, IntVect max_cells_moved,
                       bool remove_negative=true)
    {
        clearNeighbors();
        ParticleContainerType::Redistribute(lev_min, lev_max, nGrow, local,
                                            max_cells_moved, remove_negative);
    }

    void RedistributeLocal (bool remove_negative=true)
    {
        const int lev_min = 0;
        const int lev_max = 0;
        const IntVect nGrow(0);
        const IntVect local(1);
        clearNeighbors();
        this->Redistribute(lev_min, lev_max, nGrow, true, local, remove_negative);
    }

#ifdef AMREX_USE_GPU
    void fillNeighborsGPU ();
    void updateNeighborsGPU (bool boundary_neighbors_only=false);
    void clearNeighborsGPU ();
#else
    void fillNeighborsCPU (amrex::Real search_radius = amrex::Real(0));
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

    [[nodiscard]] int numRealCommComps () const { return RealCommCompStart + this->NumRealComps(); }
    [[nodiscard]] int numIntCommComps () const { return IntCommCompStart + this->NumIntComps(); }

    ///
    /// Lists of the components the same-rank ghost refresh has to copy, from
    /// the ghost flags and the slot group (see setGhostSlotGroup).
    ///
    [[nodiscard]] detail::GhostCopyPlan makeGhostCopyPlan () const
    {
        detail::GhostCopyPlan plan;
        constexpr int rcs = RealCommCompStart;
        constexpr int ics = IntCommCompStart;
        const int nrc = this->NumRealComps();
        const int nic = this->NumIntComps();
        const int nrt = nrc - T_NArrayReal;   // runtime reals
        const int nit = nic - T_NArrayInt;    // runtime ints

        // The group is used only if it lies inside the runtime components
        // with distinct key arrays, and every component of it is flagged.
        bool grp = (ghost_slot_n > 0) && (ghost_slot_nblocks > 0)
            && (ghost_slot_real0 >= 0)
            && (ghost_slot_real0 + ghost_slot_nblocks*ghost_slot_n <= nrt)
            && (ghost_slot_key_cpu >= 0) && (ghost_slot_key_id >= 0)
            && (ghost_slot_key_cpu + ghost_slot_n <= nit)
            && (ghost_slot_key_id  + ghost_slot_n <= nit)
            && (ghost_slot_key_cpu + ghost_slot_n <= ghost_slot_key_id
                || ghost_slot_key_id + ghost_slot_n <= ghost_slot_key_cpu);
        if (grp) {
            for (int k = 0; grp && k < ghost_slot_nblocks*ghost_slot_n; ++k) {
                grp = ghost_real_comp[rcs + T_NArrayReal + ghost_slot_real0 + k] != 0;
            }
            for (int k = 0; grp && k < ghost_slot_n; ++k) {
                grp = ghost_int_comp[ics + T_NArrayInt + ghost_slot_key_cpu + k] != 0
                   && ghost_int_comp[ics + T_NArrayInt + ghost_slot_key_id  + k] != 0;
            }
        }
        if (grp) {
            plan.nslots  = ghost_slot_n;
            plan.key_cpu = ghost_slot_key_cpu;
            plan.key_id  = ghost_slot_key_id;
            plan.real0   = ghost_slot_real0;
            plan.nblocks = ghost_slot_nblocks;
        }

        for (int ii = 0; ii < nrc; ++ii) {
            if (!ghost_real_comp[rcs + ii]) { continue; }
            if (ii < T_NArrayReal) {
                plan.ct_real.push_back(ii);
            } else {
                int const jr = ii - T_NArrayReal;
                if (grp && jr >= plan.real0 && jr < plan.real0 + plan.nblocks*plan.nslots) { continue; }
                plan.rt_real.push_back(jr);
            }
        }
        for (int ii = 0; ii < nic; ++ii) {
            if (!ghost_int_comp[ics + ii]) { continue; }
            if (ii < T_NArrayInt) {
                plan.ct_int.push_back(ii);
            } else {
                int const jr = ii - T_NArrayInt;
                if (grp && ((jr >= plan.key_cpu && jr < plan.key_cpu + plan.nslots) ||
                            (jr >= plan.key_id  && jr < plan.key_id  + plan.nslots))) { continue; }
                plan.rt_int.push_back(jr);
            }
        }

        // Record layout for the masked MPI path: byte offsets of the slot keys
        // and the slot of every payload component, in the order in which
        // packNeighborParticleData writes the flagged components.
        if (grp && ghost_mpi_masked && plan.nslots <= detail::MaxWireSlots) {
            std::size_t off = 0;
            plan.wire_real_off.assign(nrc, -1);
            plan.wire_int_off.assign(nic, -1);
            if constexpr (!ParticleType::is_soa_particle) {
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    if (ghost_real_comp[dir]) {
                        plan.wire_pos_off[dir] = int(off);
                        off += sizeof(typename ParticleType::RealType);
                    }
                }
                for (int ii = 0; ii < ParticleType::NReal; ++ii) {
                    if (ghost_real_comp[AMREX_SPACEDIM + ii]) { off += sizeof(typename ParticleType::RealType); }
                }
            }
            plan.wire_real_slot.assign(nrc, -1);
            for (int ii = 0; ii < nrc; ++ii) {
                if (!ghost_real_comp[rcs + ii]) { continue; }
                plan.wire_real_off[ii] = int(off);
                if (ii >= T_NArrayReal) {
                    int const jr = ii - T_NArrayReal;
                    if (jr >= plan.real0 && jr < plan.real0 + plan.nblocks*plan.nslots) {
                        plan.wire_real_slot[ii] = (jr - plan.real0) % plan.nslots;
                    }
                }
                off += sizeof(ParticleReal);
            }
            for (int ii = 0; ii < 2; ++ii) {
                if (ghost_int_comp[ii]) { plan.wire_idcpu_off[ii] = int(off); off += sizeof(int); }
            }
            if constexpr (!ParticleType::is_soa_particle) {
                for (int ii = 0; ii < ParticleType::NInt; ++ii) {
                    if (ghost_int_comp[2 + ii]) { off += sizeof(int); }
                }
            }
            plan.wire_key_cpu_off.assign(plan.nslots, -1);
            plan.wire_key_id_off.assign(plan.nslots, -1);
            for (int ii = 0; ii < nic; ++ii) {
                if (!ghost_int_comp[ics + ii]) { continue; }
                plan.wire_int_off[ii] = int(off);
                if (ii >= T_NArrayInt) {
                    int const jr = ii - T_NArrayInt;
                    if (jr >= plan.key_cpu && jr < plan.key_cpu + plan.nslots) {
                        plan.wire_key_cpu_off[jr - plan.key_cpu] = int(off);
                    }
                    if (jr >= plan.key_id && jr < plan.key_id + plan.nslots) {
                        plan.wire_key_id_off[jr - plan.key_id] = int(off);
                    }
                }
                off += sizeof(int);
            }
            plan.rec_size = int(off);
            plan.wire = true;
            if constexpr (!ParticleType::is_soa_particle) {
                // batched records: the layout above has to be the communication layout
                plan.batch = ghost_batched && (ParticleType::NReal == 0) && (ParticleType::NInt == 0)
                          && (std::size_t(off) == cdata_size);
                if (plan.batch && ghost_compact && plan.nblocks <= detail::MaxPayloadBlocks) {
                    // compact records: the same components in the same order without the
                    // payload of the slots, then the payload offset of the ghost
                    std::size_t coff = 0;
                    plan.cwire_real_off.assign(nrc, -1);
                    plan.cwire_int_off.assign(nic, -1);
                    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                        if (ghost_real_comp[dir]) {
                            plan.cwire_pos_off[dir] = int(coff);
                            coff += sizeof(typename ParticleType::RealType);
                        }
                    }
                    for (int ii = 0; ii < nrc; ++ii) {
                        if (!ghost_real_comp[rcs + ii]) { continue; }
                        if (plan.wire_real_slot[ii] >= 0) { continue; }
                        plan.cwire_real_off[ii] = int(coff);
                        coff += sizeof(ParticleReal);
                    }
                    for (int ii = 0; ii < 2; ++ii) {
                        if (ghost_int_comp[ii]) { plan.cwire_idcpu_off[ii] = int(coff); coff += sizeof(int); }
                    }
                    plan.cwire_key_cpu_off.assign(plan.nslots, -1);
                    plan.cwire_key_id_off.assign(plan.nslots, -1);
                    for (int ii = 0; ii < nic; ++ii) {
                        if (!ghost_int_comp[ics + ii]) { continue; }
                        plan.cwire_int_off[ii] = int(coff);
                        if (ii >= T_NArrayInt) {
                            int const jr = ii - T_NArrayInt;
                            if (jr >= plan.key_cpu && jr < plan.key_cpu + plan.nslots) {
                                plan.cwire_key_cpu_off[jr - plan.key_cpu] = int(coff);
                            }
                            if (jr >= plan.key_id && jr < plan.key_id + plan.nslots) {
                                plan.cwire_key_id_off[jr - plan.key_id] = int(coff);
                            }
                        }
                        coff += sizeof(int);
                    }
                    plan.cpoff_off  = int(coff);
                    coff += sizeof(int);
                    plan.chead_size = int(coff);
                    plan.cpay_bytes = plan.nblocks * int(sizeof(ParticleReal));
                    plan.compact = true;
                }
            }
        }
        return plan;
    }

    void cacheNeighborInfo (amrex::Real search_radius = amrex::Real(0));

    ///
    /// This builds the internal mask data structure used for looking up neighbors.
    /// When search_radius > 0, the mask grow is sized to cover that physical
    /// distance per level; otherwise the constructor's m_num_neighbor_cells
    /// is used.
    ///
    void BuildMasks (amrex::Real search_radius = amrex::Real(0));

    ///
    /// Are the masks computed by the above function still valid?
    /// When search_radius > 0, masks are additionally invalid if their
    /// per-level grow is smaller than ceil(search_radius / min(dx_lev)).
    ///
    bool areMasksValid (amrex::Real search_radius = amrex::Real(0));

    void GetNeighborCommTags (amrex::Real search_radius = amrex::Real(0));

    void GetCommTagsBox (Vector<NeighborCommTag>& tags, int lev, const Box& in_box,
                         amrex::Real search_radius = amrex::Real(0));

    void resizeContainers (int num_levels);

    void initializeCommComps ();

    void calcCommSize ();

    ///
    /// Perform the MPI communication necessary to fill neighbor buffers
    ///
    void fillNeighborsMPI (bool reuse_rcv_counts, detail::GhostCopyPlan const& plan,
                           bool direct = false, bool batched = false);

    void sumNeighborsMPI (std::map<int, Vector<char> >& not_ours,
                          int real_start_comp, int real_num_comp,
                          int int_start_comp, int int_num_comp);

    ///
    /// Perform handshake to figure out how many bytes each proc should receive
    ///
    void getRcvCountsMPI ();

    template <typename P>
    void getNeighborTags (Vector<NeighborCopyTag>& tags,
                          const P& p,
                          int nGrow, const NeighborCopyTag& src_tag, const MyParIter& pti,
                          amrex::Real search_radius = amrex::Real(0));

    template <typename P>
    void getNeighborTags (Vector<NeighborCopyTag>& tags,
                          const P& p,
                          const IntVect& nGrow, const NeighborCopyTag& src_tag, const MyParIter& pti,
                          amrex::Real search_radius = amrex::Real(0));

    IntVect computeRefFac (int src_lev, int lev);

    Vector<std::map<PairIndex, Vector<InverseCopyTag> > > inverse_tags;
    Vector<std::map<PairIndex, ParticleTile> > neighbors;
    Vector<std::map<PairIndex, IntVector> >      neighbor_list;
    static constexpr size_t pdata_size = sizeof(ParticleType);

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

    // Slot group and mode of the same-rank ghost refresh (CPU)
    int ghost_slot_n       = 0;
    int ghost_slot_key_cpu = 0;
    int ghost_slot_key_id  = 0;
    int ghost_slot_real0   = 0;
    int ghost_slot_nblocks = 0;
    bool ghost_copy_masked = false;
    bool ghost_mpi_masked  = false;

    // Options of the ghost refresh on the MPI path (CPU): persistent receive buffer and plan,
    // ghosts written straight into the tiles, records packed per component
    bool ghost_persistent_recv = false;
    bool ghost_direct          = false;
    bool ghost_batched         = false;
    bool ghost_isend           = false;
    // Compact records (setGhostCompact): switch, whether the layout of the last rebuild is the compact one,
    // bytes of the head part of the message to each rank (the payload area follows), and the state of the
    // payload area of each rank in the current update
    bool ghost_compact         = false;
    bool ghost_compact_layout  = false;
    std::map<int, std::size_t> ghost_head_bytes;
    Vector<detail::GhostSendInfo> ghost_send_info;
    bool ghost_plan_valid      = false;
    bool ghost_switches_told   = false;
    Vector<char> ghost_recv_buf;
    detail::GhostCopyPlan ghost_plan_cache;
    detail::GhostBatchScratch ghost_scratch;

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
            if (grid_id != other.grid_id) { return grid_id < other.grid_id; }
            if (box     != other.box    ) { return box < other.box; }
            AMREX_D_TERM(
            if (periodic_shift[0] != other.periodic_shift[0])
                { return periodic_shift[0] < other.periodic_shift[0]; },
            if (periodic_shift[1] != other.periodic_shift[1])
                { return periodic_shift[1] < other.periodic_shift[1]; },
            if (periodic_shift[2] != other.periodic_shift[2])
                { return periodic_shift[2] < other.periodic_shift[2]; }
                )
            return false;
        }
    };

    // These are used to keep track of which particles need to be ghosted to which grids
    bool m_neighbor_mask_initialized = false;
    amrex::BoxArray m_cached_neighbor_ba;
    amrex::DistributionMapping m_cached_neighbor_dm;

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
    amrex::PODVector<char, PolymorphicArenaAllocator<char> > rcv_buffer;

    Gpu::PinnedVector<char> pinned_snd_buffer;
    Gpu::PinnedVector<char> pinned_rcv_buffer;
#endif

    NeighborListContainerType m_neighbor_list;

    Vector<std::map<std::pair<int, int>, amrex::Gpu::DeviceVector<int> > > m_boundary_particle_ids;

    [[nodiscard]] bool hasNeighbors() const { return m_has_neighbors; }

    bool m_has_neighbors = false;
};

template <int T_NStructReal, int T_NStructInt, int T_NArrayReal=0, int T_NArrayInt=0>
using NeighborParticleContainer =
    NeighborParticleContainer_impl<Particle<T_NStructReal, T_NStructInt>, T_NArrayReal, T_NArrayInt>;

template <int T_NArrayReal, int T_NArrayInt,
          template<class> class Allocator=DefaultAllocator,
          class CellAssignor=DefaultAssignor>
using NeighborParticleContainerPureSoA =
    NeighborParticleContainer_impl<SoAParticle<T_NArrayReal, T_NArrayInt>,
                                   T_NArrayReal, T_NArrayInt, Allocator, CellAssignor>;

}

#include "AMReX_NeighborParticlesI.H"

#ifdef AMREX_USE_GPU
#include "AMReX_NeighborParticlesGPUImpl.H"
#else
#include "AMReX_NeighborParticlesCPUImpl.H"
#endif

#endif // _NEIGHBORPARTICLES_H_
