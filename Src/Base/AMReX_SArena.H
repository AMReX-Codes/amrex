#ifndef AMREX_SARENA_H_
#define AMREX_SARENA_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_SArena.H
 * \brief Stream-ordered memory arena that defers frees until after pending GPU kernels complete.
 */

#include <AMReX_Arena.H>


namespace amrex {
/**
 * \class amrex::SArena
 * \ingroup amrex_memory
 * \brief A STREAM-ordered memory arena.
*
* This Arena is implemented as a wrapper around amrex::The_Arena()
* where the free() function calls amrex::Gpu::freeAsync().
* This allows memory to be used by GPU kernels that were launched before free() was called
* but may execute afterward, without calling amrex::Gpu::streamSynchronize().
* These kernels need to use the same GPU stream that is active when free() is called,
* otherwise explicit synchronization is needed.
*
 * This is achieved by holding freed memory in a temporary buffer until the next time the
 * GPU stream is synchronized. In case too much memory is held up in this buffer or if The_Arena
 * needs more memory, the stream is automatically synchronized to clear the buffer.
 */
class SArena
    :
    public Arena
{
public:
    /**
    * \brief Allocates dynamic memory from the arena of size sz.
    * Returns a pointer to this memory.
    * The memory can be used immediately by both the CPU and GPU with any stream.
    * This function may synchronize all GPU streams to free up memory.
    */
    [[nodiscard]] void* alloc (std::size_t sz) final;
    /**
    * \brief Deletes the arena pointed to by pt.
    * After this function is called, the memory can still be used by already launched
    * GPU kernels on the currently active stream.
    * This function may synchronize the currently active GPU stream to free up memory.
    */
    void free (void* pt) final;

    //! Return true; SArena allocations are accessible from the device.
    [[nodiscard]] bool isDeviceAccessible () const final;
    //! Return true; SArena allocations reside in host memory.
    [[nodiscard]] bool isHostAccessible () const final;

    //! Return false; SArena does not allocate managed (unified) memory.
    [[nodiscard]] bool isManaged () const final;
    //! Return false because SArena is backed by host memory.
    [[nodiscard]] bool isDevice () const final;
    //! Return false; SArena allocations are not pinned host buffers.
    [[nodiscard]] bool isPinned () const final;

#ifdef AMREX_USE_GPU
    //! Is this a GPU stream ordered memory allocator?
    [[nodiscard]] bool isStreamOrderedArena () const final { return true; }
#endif
};

}

#endif
