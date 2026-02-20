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
    //! Construct with a memory-pool release threshold (bytes).
    PArena (Long release_threshold);
    PArena (const PArena& rhs) = delete;
    PArena (PArena&& rhs) = delete;
    PArena& operator= (const PArena& rhs) = delete;
    PArena& operator= (PArena&& rhs) = delete;
    ~PArena () override;

    [[nodiscard]] void* alloc (std::size_t nbytes) final;
    void free (void* p) final;

    [[nodiscard]] bool isDeviceAccessible () const final;
    [[nodiscard]] bool isHostAccessible () const final;

    [[nodiscard]] bool isManaged () const final;
    [[nodiscard]] bool isDevice () const final;
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
