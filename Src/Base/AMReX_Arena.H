#ifndef AMREX_ARENA_H_
#define AMREX_ARENA_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_Arena.H
 * \brief Memory arena base class and global arena accessors.
 */

#include <AMReX_BLassert.H>
#include <AMReX_INT.H>
#ifdef AMREX_USE_GPU
#include <AMReX_GpuControl.H>
#endif

#ifdef AMREX_TINY_PROFILING
#include <AMReX_TinyProfiler.H>
#else
namespace amrex {
    struct MemStat {};
}
#endif

#include <cstddef>
#include <cstdlib>
#include <limits>
#include <map>
#include <mutex>
#include <unordered_map>
#include <utility>

namespace amrex {

/**
 * \ingroup amrex_memory
 * \brief Return the smallest multiple of \p align_requirement that is >= \p size.
 */
inline std::size_t aligned_size (std::size_t align_requirement, std::size_t size) noexcept
{
    return ((size + (align_requirement-1)) / align_requirement) * align_requirement;
}

/// \ingroup amrex_memory
//! Return whether the address \p p is aligned to \p alignment bytes.
inline bool is_aligned (const void* p, std::size_t alignment) noexcept
{
    auto* q = const_cast<void*>(p);
    auto space = alignment;
    return std::align(alignment, alignment, q, space) == p;
}

class Arena;

//! \ingroup amrex_memory
//! Return a pointer to `The_Arena`, the main memory arena.
Arena* The_Arena ();
//! \ingroup amrex_memory
//! Return a pointer to `The_Async_Arena`, the async-safe GPU device arena.
Arena* The_Async_Arena ();
//! \ingroup amrex_memory
//! Return a pointer to `The_Device_Arena`, the GPU device memory arena. It
//! may be an alias to The_Arena.
Arena* The_Device_Arena ();
//! \ingroup amrex_memory
//! Return a pointer to `The_Managed_Arena`, the GPU managed memory
//! arena. It may be an alias to The_Arena.
Arena* The_Managed_Arena ();
//! \ingroup amrex_memory
//! Return a pointer to `The_Pinned_Arena`, the pinned memory arena.
Arena* The_Pinned_Arena ();
//! \ingroup amrex_memory
//! Return a pointer to `The_Comms_Arena`, the memory arena used by
//! communication buffers.
Arena* The_Comms_Arena ();
//! \ingroup amrex_memory
//! Return a pointer to `The_Cpu_Arena`, the CPU memory arena.
Arena* The_Cpu_Arena ();

//! Return a process-lifetime no-op arena used once other arenas are torn down.
//! Note: still unsafe as a captured deallocator in foreign static objects,
//! because C++ teardown order can destroy the null arena before they run.
Arena* The_Null_Arena ();

/**
 * \ingroup amrex_memory
 * \brief Configuration options for constructing an Arena.
 *
 * Uses a builder pattern: each setter returns \c *this to allow chaining.
 */
struct ArenaInfo
{
    Long release_threshold = std::numeric_limits<Long>::max(); ///< When usage would exceed this, cached unused blocks are released first.
    bool use_cpu_memory = false;              ///< Allocate CPU (host) memory instead of GPU memory.
    bool device_use_managed_memory = true;    ///< Use CUDA/HIP managed (unified) memory for GPU allocations.
    bool device_set_readonly = false;         ///< Advise the runtime that the managed memory is read-only on the device.
    bool device_set_preferred = false;        ///< Advise the runtime to prefer the device for managed memory.
    bool device_use_hostalloc = false;        ///< Use host-pinned memory (cudaMallocHost/hipHostMalloc).
    bool defragmentation = false;             ///< Enable defragmentation in the arena.
    //! Set the usage threshold that triggers release of cached unused blocks before growing.
    ArenaInfo& SetReleaseThreshold (Long rt) noexcept {
        release_threshold = rt;
        return *this;
    }
    //! Configure the arena to use non-managed GPU device memory.
    ArenaInfo& SetDeviceMemory () noexcept {
        device_use_managed_memory = false;
        device_use_hostalloc = false;
        return *this;
    }
    //! Advise the runtime that the managed memory will be read-only on the device.
    ArenaInfo& SetReadOnly () noexcept {
        BL_ASSERT(device_use_managed_memory);
        device_set_readonly = true;
        return *this;
    }
    //! Advise the runtime to prefer device residence for the managed memory.
    ArenaInfo& SetPreferred () noexcept {
        BL_ASSERT(device_use_managed_memory);
        device_set_preferred = true;
        return *this;
    }
    //! Configure the arena to use host-pinned memory.
    ArenaInfo& SetHostAlloc () noexcept {
        device_use_hostalloc = true;
        device_use_managed_memory = false;
        return *this;
    }
    //! Configure the arena to use CPU (host) memory.
    ArenaInfo& SetCpuMemory () noexcept {
        use_cpu_memory = true;
        device_use_managed_memory = false;
        device_set_readonly = false;
        device_set_preferred = false;
        device_use_hostalloc = false;
        return *this;
    }
    //! Enable or disable defragmentation in the arena.
    ArenaInfo& SetDefragmentation (bool b) noexcept {
        defragmentation = b;
        return *this;
    }
};

/**
 * \class amrex::Arena
 * \ingroup amrex_memory
 * \brief Abstract base class for memory arenas.
 *
 * Subclasses manage a pool of memory (CPU, GPU device, managed, or pinned)
 * and expose alloc()/free() for raw byte allocation and deallocation.
 */
class Arena
{
public:

    virtual ~Arena () = default;
    Arena () noexcept = default;

    Arena (const Arena& rhs) = delete;
    Arena (Arena&& rhs) = delete;
    Arena& operator= (const Arena& rhs) = delete;
    Arena& operator= (Arena&& rhs) = delete;

    /**
     * \brief Allocate \p sz bytes from this arena.
     * \param sz Number of bytes to allocate.
     * \return Pointer to the allocated memory.
     */
    [[nodiscard]] virtual void* alloc (std::size_t sz) = 0;

    /**
     * \brief Allocate memory in-place if possible.
     *
     * In this implementation of the virtual function, a new allocation is
     * always performed, and it returns the new pointer and allocation size.
     * Note that this function does NOT free the existing allocation. So
     * it's the caller's responsibility to free the original pointer when
     * appropriate.
     *
     * \param pt Existing allocation to be reused when possible.
     * \param szmin Minimum acceptable size for reusing \p pt in place.
     * \param szmax Requested size for the returned allocation.
     * \return Pointer/size pair describing the resulting allocation.
     */
    [[nodiscard]] virtual std::pair<void*,std::size_t>
    alloc_in_place (void* pt, std::size_t szmin, std::size_t szmax)
    {
        ignore_unused(pt, szmin);
        auto* p = alloc(szmax);
        return std::make_pair(p, szmax);
    }

    /**
     * \brief Shrink allocation size in-place.
     *
     * In this implementation of the virtual function, a new allocation is
     * always performed, and it returns the new pointer. Note that this
     * function does NOT free the existing allocation, so it's the caller's
     * responsibility to free the original pointer when appropriate.
     *
     * \param pt Existing allocation to be shrunk.
     * \param sz Target size for the new allocation.
     * \return Pointer to the resized allocation.
     */
    [[nodiscard]] virtual void*
    shrink_in_place (void* pt, std::size_t sz)
    {
        ignore_unused(pt);
        return alloc(sz);
    }

    //! Free a previously allocated block pointed to by \p pt.
    virtual void free (void* pt) = 0;

    /**
     * \brief Release unused cached memory back to the system.
     * \return Number of bytes freed.
     */
    virtual std::size_t freeUnused () { return 0; }

    //! Check whether the memory is accessible on the device. Note that
    //! isDeviceAccessible and isHostAccessible can both be true.
    [[nodiscard]] virtual bool isDeviceAccessible () const;

    //! Check whether the memory is accessible on the host. Note that
    //! isDeviceAccessible and isHostAccessible can both be true.
    [[nodiscard]] virtual bool isHostAccessible () const;

    /**
     * \brief Check whether it is managed GPU memory
     *
     * Note that isManaged, isDevice and isPinned are mutually exclusive.
     * For memory allocated by cudaMalloc* etc. in GPU builds, one of them
     * returns true. For CPU builds, they are always false.
     */
    [[nodiscard]] virtual bool isManaged () const;

    /**
     * \brief Check whether it is non-managed GPU device memory
     *
     * Note that isManaged, isDevice and isPinned are mutually exclusive.
     * For memory allocated by cudaMalloc* etc. in GPU builds, one of them
     * returns true. For CPU builds, they are always false.
     */
    [[nodiscard]] virtual bool isDevice () const;

    /**
     * \brief Check whether it is pinned host memory
     *
     * Note that isManaged, isDevice and isPinned are mutually exclusive.
     * For memory allocated by cudaMalloc* etc. in GPU builds, one of them
     * returns true. For CPU builds, they are always false.
     */
    [[nodiscard]] virtual bool isPinned () const;

    /**
     * \brief Does the device have enough free memory for allocating this
     * much memory?  For CPU builds, this always return true.  This is not a
     * const function because it may attempt to release memory back to the
     * system.
     * \param sz Requested allocation size in bytes.
     */
    [[nodiscard]] virtual bool hasFreeDeviceMemory (std::size_t sz);

    /**
     * \brief Add this Arena to the list of Arenas that are profiled by TinyProfiler.
     * \param memory_name The name of this arena in the TinyProfiler output.
     */
    void registerForProfiling (const std::string& memory_name);

    /**
     * \brief Remove this Arena from the list of Arenas that are profiled by TinyProfiler.
     * This is equivalent to destructing and re-constructing the Arena
     */
    void deregisterFromProfiling ();

#ifdef AMREX_USE_GPU
    //! Is this GPU stream ordered memory allocator?
    [[nodiscard]] virtual bool isStreamOrderedArena () const { return false; }
    virtual void streamOrderedFree (void* pt, gpuStream_t stream);
#endif

    /**
     * \brief Return the smallest multiple of Arena::align_size that is >= \p sz bytes.
     */
    static std::size_t align (std::size_t sz);

    /**
     * \brief Initialize the Arena subsystem (used internally by AMReX).
     * \param minimal When true, perform only the minimal setup required.
     */
    static void Initialize (bool minimal);

    /**
     * \brief Print memory usage information of all arenas
     *
     * This function reports memory usage statistics for all `CArena`-based
     * AMReX arenas. It contains MPI collective calls and therefore must be
     * called by all processes. The output is printed on the I/O process
     * only.
     *
     * Depending on the value of `print_max_usage`, the function prints
     * either the current memory usage (when false) or the maximum memory
     * usage recorded (when true).
     *
     * \param print_max_usage whether to print the maximum usage or the
     *        current usage
     */
    static void PrintUsage (bool print_max_usage=false);

    /**
     * \brief Print memory usage information of all arenas to a given output
     * stream
     *
     * This function reports memory usage statistics for all `CArena`-based
     * AMReX arenas. Unlike `PrintUsage`, this does not contain MPI calls
     * and therefore may be called independently on individual processes.
     *
     * \param os output stream
     * \param space indentation prepended to printed lines
     */
    static void PrintUsageToStream (std::ostream& os, std::string const& space);

    /**
     * \brief Print memory usage information of all arenas to given file
     *
     * This function reports memory usage statistics for all `CArena`-based
     * AMReX arenas. The output is appended to a file named
     * `{filename}.{rank}`, where `rank` is the MPI rank of the process.
     * Unlike `PrintUsage`, this does not contain MPI calls and therefore
     * may be called independently on individual processes.
     *
     * \param filename base name of the output file
     * \param message message printed at the beginning
     */
    static void PrintUsageToFiles (std::string const& filename, std::string const& message);

    static void Finalize (); ///< Used internally by amrex

    /// The alignment of allocated memory
    /// 256 bytes is the default alignment that cudaMalloc uses.
    /// If the native allocation function were to give memory that
    /// has a smaller alignment than align_size, memory from
    /// the Arena will also only have this smaller alignment.
    /// In practice, native allocation functions have >= 256 bytes
    /// of alignment. In all cases, memory should be at least
    /// aligned to std::max_align_t, usually 16.
    static const std::size_t align_size = 256;

    /**
     *  \brief Return the ArenaInfo object for querying
     */
    [[nodiscard]] const ArenaInfo& arenaInfo () const { return arena_info; }

    /// Reset the maximum usage counter. This is no-op in this
    /// implementation of the virtual function.
    virtual void ResetMaxUsageCounter () {}

    /** Whether the Arena subsystem is initialized. */
    [[nodiscard]] static bool IsInitialized ();

protected:

    ArenaInfo arena_info;

    virtual std::size_t freeUnused_protected () { return 0; }
    void* allocate_system (std::size_t nbytes);
    void deallocate_system (void* p, std::size_t nbytes);
    static void out_of_memory_abort (std::string const& memory_type, std::size_t nbytes,
                                     std::string const& error_msg);

    /// \cond DOXYGEN_IGNORE
    struct ArenaProfiler {
        //! If this arena is profiled by TinyProfiler
        bool m_do_profiling = false;
        //! Mutex for the profiling
        std::mutex m_arena_profiler_mutex;
        //! Data structure used for profiling with TinyProfiler
        std::map<std::string, MemStat> m_profiling_stats;
        //! Track the currently allocated memory, not used by CArena
        std::unordered_map<void*, std::pair<MemStat*, std::size_t>> m_currently_allocated;

        ~ArenaProfiler ();
        ArenaProfiler () noexcept = default;
        ArenaProfiler (const ArenaProfiler& rhs) = delete;
        ArenaProfiler (ArenaProfiler&& rhs) = delete;
        ArenaProfiler& operator= (const ArenaProfiler& rhs) = delete;
        ArenaProfiler& operator= (ArenaProfiler&& rhs) = delete;

        void profile_alloc (void* ptr, std::size_t nbytes);

        void profile_free (void* ptr);

    } m_profiler;
    /// \endcond
};

}

#endif /*BL_ARENA_H*/
