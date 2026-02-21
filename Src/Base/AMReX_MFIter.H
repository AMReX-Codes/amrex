#ifndef BL_MFITER_H_
#define BL_MFITER_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_MFIter.H
 * \brief Iterator for looping over tiles and boxes of FabArray-based containers.
 */

#include <AMReX_FabArrayBase.H>

#include <memory>

namespace amrex {

/// \ingroup amrex_fab_containers
//! Return \c true when tiling should be used: always on CPU, only outside GPU launch regions on GPU.
#ifdef AMREX_USE_GPU
    inline bool TilingIfNotGPU () noexcept { return Gpu::notInLaunchRegion(); }
#else
    constexpr bool TilingIfNotGPU () noexcept { return true; }
#endif

template<class T> class FabArray;

/**
 * \struct amrex::MFItInfo
 * \ingroup amrex_fab_containers
 * \brief Configuration options passed to MFIter constructors.
 *
 * Controls tiling, dynamic scheduling, GPU stream count, and GPU device
 * synchronization behaviour.  All setters return \c *this for chaining.
 */
struct MFItInfo
{
    bool do_tiling{false};      //!< Enable tiling when \c true.
    bool dynamic{false};        //!< Use dynamic (work-stealing) scheduling when \c true.
    bool device_sync = true;    //!< Synchronize device before and after iteration when \c true.
    bool device_sync_pre = true;   //!< Synchronize device before iteration starts.
    bool device_sync_post = true;  //!< Synchronize device after iteration ends.
    int  num_streams = 1;       //!< Number of GPU streams to use.
    IntVect tilesize;           //!< Tile size (zero means use default).
    //! Default-construct with GPU-aware defaults.
    MFItInfo () noexcept
        :  device_sync(!Gpu::inNoSyncRegion()), num_streams(Gpu::numGpuStreams()),
          tilesize(IntVect::TheZeroVector()) {}
    //! Enable tiling with tile size \p ts (defaults to \c FabArrayBase::mfiter_tile_size).
    MFItInfo& EnableTiling (const IntVect& ts = FabArrayBase::mfiter_tile_size) noexcept {
        do_tiling = true;
        tilesize = ts;
        return *this;
    }
    //! Set tiling on or off.
    MFItInfo& SetTiling (bool f) noexcept {
        do_tiling = f;
        return *this;
    }
    //! Enable or disable dynamic (work-stealing) scheduling.
    MFItInfo& SetDynamic (bool f) noexcept {
        dynamic = f;
        return *this;
    }
    //! Disable all GPU device synchronization around the iteration.
    MFItInfo& DisableDeviceSync () noexcept {
        device_sync = false;
        return *this;
    }
    //! Enable or disable all GPU device synchronization around the iteration.
    MFItInfo& SetDeviceSync (bool f) noexcept {
        device_sync = f;
        return *this;
    }
    //! Disable the GPU device synchronization performed before iteration starts.
    MFItInfo& DisableDeviceSyncPre () noexcept {
        device_sync_pre = false;
        return *this;
    }
    //! Enable or disable the GPU device synchronization before iteration starts.
    MFItInfo& SetDeviceSyncPre (bool f) noexcept {
        device_sync_pre = f;
        return *this;
    }
    //! Disable the GPU device synchronization performed after iteration ends.
    MFItInfo& DisableDeviceSyncPost () noexcept {
        device_sync_post = false;
        return *this;
    }
    //! Enable or disable the GPU device synchronization after iteration ends.
    MFItInfo& SetDeviceSyncPost (bool f) noexcept {
        device_sync_post = f;
        return *this;
    }
    //! Set the number of GPU streams to use during iteration.
    MFItInfo& SetNumStreams (int n) noexcept {
        num_streams = n;
        return *this;
    }
    //! Restrict iteration to the default GPU stream (sets \c num_streams to 1).
    MFItInfo& UseDefaultStream () noexcept {
        num_streams = 1;
        return *this;
    }
};

/**
 * \class amrex::MFIter
 * \ingroup amrex_multifab_containers
 * \brief Iterator for looping over tiles and boxes of \ref amrex::FabArray
 * based containers.
 */
class MFIter
{
public:
    //! Flags that specify tiling modes. All these flags are off by default.
    enum Flags {
        //!Tiling: Enabling the tiling mode
        Tiling        = 0x01,
        /**
        * \brief AllBoxes: If on, all threads/workers loop over all boxes without tiling.
        * This essentially loops over indexMap.
        * Note that many functions won't work with this.
        */
        AllBoxes      = 0x02,
        //! NoTeamBarrier: This option is for Team only. If on, there is no barrier in MFIter dtor.
        NoTeamBarrier = 0x04
    };

    /**
    * \brief The default constructor does not enable tiling (flags_ is set to 0 by default).
    * However, the tiling mode can be enabled by explicitly turning the tiling flag on (i.e. flags_ = Tiling).
    * The tile size in this case is defined by FabArrayBase::mfiter_tile_size
    */
    explicit MFIter (const FabArrayBase& fabarray,
                     unsigned char       flags_=0);

    //! Enable tiling with the default tile size, which is defined by FabArrayBase::mfiter_tile_size
    MFIter (const FabArrayBase& fabarray,
            bool                do_tiling);

    //! Enable tiling with explicit tile size and flags (See type Flags for more information)
    MFIter (const FabArrayBase& fabarray,
            const IntVect&      tilesize,
            unsigned char       flags_=0);

    //! Construct over a BoxArray/DistributionMapping pair with optional flags.
    MFIter (const BoxArray& ba, const DistributionMapping& dm, unsigned char flags_=0);

    //! Construct over a BoxArray/DistributionMapping pair with tiling on or off.
    MFIter (const BoxArray& ba, const DistributionMapping& dm, bool do_tiling);

    //! Construct over a BoxArray/DistributionMapping pair with explicit tile size and flags.
    MFIter (const BoxArray& ba, const DistributionMapping& dm,
            const IntVect& tilesize, unsigned char flags_=0);

    //! Construct over a FabArrayBase with the full \c MFItInfo configuration.
    MFIter (const FabArrayBase& fabarray, const MFItInfo& info);

    //! Construct over a BoxArray/DistributionMapping pair with the full \c MFItInfo configuration.
    MFIter (const BoxArray& ba, const DistributionMapping& dm, const MFItInfo& info);

    MFIter (MFIter&& rhs) = default;
    MFIter (MFIter const&) = delete;
    MFIter& operator= (MFIter const&) = delete;
    MFIter& operator= (MFIter &&) = delete;

    // dtor
    ~MFIter ();

    //! Return the tile Box at the current index.
    [[nodiscard]] Box tilebox () const noexcept;

    //! Return the tilebox with provided nodal flag
    [[nodiscard]] Box tilebox (const IntVect& nodal) const noexcept;

    //! Return the tilebox with provided nodal flag and grown cells
    [[nodiscard]] Box tilebox (const IntVect& nodal, const IntVect& ngrow) const noexcept;

    //! Return the dir-nodal (or all nodal if dir<0) Box at the current index.
    [[nodiscard]] Box nodaltilebox (int dir=-1) const noexcept;

    //! Return the tile box at the current index grown to include ghost cells.
    [[nodiscard]] Box growntilebox (int ng=-1000000) const noexcept;

    [[nodiscard]] Box growntilebox (const IntVect& ng) const noexcept;

    //! Return the dir-nodal (or all nodal if dir<0) box grown to include ghost cells.
    [[nodiscard]] Box grownnodaltilebox (int dir=-1, int ng=-1000000) const noexcept;

    [[nodiscard]] Box grownnodaltilebox (int dir, const IntVect& ng) const noexcept;

    //! Return the valid Box in which the current tile resides.
    [[nodiscard]] Box validbox () const noexcept { return fabArray->box((*index_map)[currentIndex]); }

    //! Return the Box of the FAB at which we currently point.
    [[nodiscard]] Box fabbox () const noexcept { return fabArray->fabbox((*index_map)[currentIndex]); }

    //! Increment iterator to the next tile we own.
    void operator++ () noexcept;

    //! Is the iterator valid i.e. is it associated with a FAB?
    [[nodiscard]] bool isValid () const noexcept { return currentIndex < endIndex; }

    //! The index into the underlying BoxArray of the current FAB.
    [[nodiscard]] int index () const noexcept { return (*index_map)[currentIndex]; }

    //! The number of indices.
    [[nodiscard]] int length () const noexcept { return (endIndex - beginIndex); }

    //! The current local tile index in the current grid;
    [[nodiscard]] int LocalTileIndex () const noexcept {return local_tile_index_map ? (*local_tile_index_map)[currentIndex] : 0;}

    //! The the number of tiles in the current grid;
    [[nodiscard]] int numLocalTiles() const noexcept {return num_local_tiles ? (*num_local_tiles)[currentIndex] : 1;}

    /**
    * \brief Return local index into the vector of fab pointers, m_fabs_v
    * When AllBoxes is on, local_index_map is a nullptr and local index is current index.
    */
    [[nodiscard]] int LocalIndex () const noexcept { return local_index_map ? (*local_index_map)[currentIndex] : currentIndex; }

    //! Constant reference to FabArray over which we're iterating.
    [[nodiscard]] const FabArrayBase& theFabArrayBase () const noexcept { return *fabArray; }

    //! Return the internal sequential tile index (for advanced use).
    [[nodiscard]] int tileIndex () const noexcept {return currentIndex;}

    //! Return the DistributionMapping of the underlying FabArray.
    [[nodiscard]] const DistributionMapping& DistributionMap () const noexcept { return fabArray->DistributionMap(); }

    //! Allow or disallow nested MFIters; returns the old value.
    static int allowMultipleMFIters (int allow);

    //! Return the current MFIter nesting depth.
    static int currentDepth ();

    //! Finalize the iterator early (otherwise called automatically in the destructor).
    void Finalize ();

protected:

    std::unique_ptr<FabArrayBase> m_fa;  //!< This must be the first member!

    const FabArrayBase* fabArray;

    IntVect tile_size;

    unsigned char flags;
    int           currentIndex;
    int           beginIndex;
    int           endIndex;
    int           streams;
    IndexType     typ;

    bool          dynamic;
    bool          finalized = false;

    struct DeviceSync {
        DeviceSync (bool f) : flag(f) {}
        DeviceSync (DeviceSync&& rhs)  noexcept : flag(std::exchange(rhs.flag,false)) {}
        ~DeviceSync () = default;
        DeviceSync (DeviceSync const&) = delete;
        DeviceSync& operator= (DeviceSync const&) = delete;
        DeviceSync& operator= (DeviceSync &&) = delete;
        explicit operator bool() const noexcept { return flag; }
        bool flag = true;
    };
    DeviceSync device_sync;
    DeviceSync device_sync_pre{true};
    DeviceSync device_sync_post{true};

    const Vector<int>* index_map;
    const Vector<int>* local_index_map;
    const Vector<Box>* tile_array;
    const Vector<int>* local_tile_index_map;
    const Vector<int>* num_local_tiles;

    static AMREX_EXPORT int nextDynamicIndex;
    static AMREX_EXPORT int depth;
    static AMREX_EXPORT int allow_multiple_mfiters;

    void Initialize ();
};

//! Is it safe to have these two MultiFabs in the same MFiter?
//! True means safe; false means maybe.
inline bool isMFIterSafe (const FabArrayBase& x, const FabArrayBase& y) {
    return x.DistributionMap() == y.DistributionMap()
        && BoxArray::SameRefs(x.boxArray(), y.boxArray());
}

}

#endif
