#ifndef AMREX_PARENA_H_
#define AMREX_PARENA_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_PArena.H
 * \brief CUDA stream-ordered memory arena backed by cudaMemPool / hipMemPool.
 */

#include <AMReX_Arena.H>
#include <AMReX_GpuControl.H>

#ifdef AMREX_USE_CUDA
#include <cuda.h>
#endif

namespace amrex {

/**
 * \class amrex::PArena
 * \ingroup amrex_memory
 * \brief Stream-ordered memory arena backed by \c cudaMemPool / \c hipMemPool.
 *
 * Uses the GPU stream-ordered memory allocator when available; falls back to
 * \c The_Arena() otherwise.
 */
class PArena
    :
    public Arena
{
public:
    /**
     * \brief Construct with a memory-pool release threshold (bytes).
     *
     * \param release_threshold Bytes of cached allocations the pool may keep before returning
     *                          memory to the driver/runtime.
     */
    PArena (Long release_threshold);
    PArena (const PArena& rhs) = delete;
    PArena (PArena&& rhs) = delete;
    PArena& operator= (const PArena& rhs) = delete;
    PArena& operator= (PArena&& rhs) = delete;
    //! Destroy the arena and release any outstanding pool reservations.
    ~PArena () override;

    /**
     * \brief Allocate \p nbytes from the underlying GPU memory pool.
     *
     * May synchronize the current stream if the pool needs to reclaim memory.
     *
     * \param nbytes Number of bytes to allocate.
     * \return Pointer to the allocated buffer.
     */
    [[nodiscard]] void* alloc (std::size_t nbytes) final;

    /**
     * \brief Return the allocation at \p p to the memory pool.
     *
     * The pointer must have been obtained from alloc().
     *
     * \param p Pointer to the allocation to free.
     */
    void free (void* p) final;

    //! Return true because PArena allocations are accessible from the device.
    [[nodiscard]] bool isDeviceAccessible () const final;
    //! Return false; PArena allocations reside in device memory.
    [[nodiscard]] bool isHostAccessible () const final;

    //! Return false; PArena does not provide managed memory.
    [[nodiscard]] bool isManaged () const final;
    //! Return true; PArena allocates from device memory pools.
    [[nodiscard]] bool isDevice () const final;
    //! Return false; PArena does not allocate pinned host memory.
    [[nodiscard]] bool isPinned () const final;

#ifdef AMREX_USE_GPU
    //! Is this CUDA stream ordered memory allocator?
    [[nodiscard]] bool isStreamOrderedArena () const final { return true; }
#endif

#ifdef AMREX_GPU_STREAM_ALLOC_SUPPORT
private:
    AMREX_HIP_OR_CUDA(
        hipMemPool_t m_pool;,
        cudaMemPool_t m_pool;
    )
    AMREX_HIP_OR_CUDA(
        uint64_t m_old_release_threshold;,
        cuuint64_t m_old_release_threshold;
    )
#endif
};

}
#endif
