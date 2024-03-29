#ifndef AMREX_PARENA_H_
#define AMREX_PARENA_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>
#include <AMReX_GpuControl.H>

#ifdef AMREX_USE_CUDA
#include <cuda.h>
#endif

namespace amrex {

/**
* \brief This arena uses CUDA stream-ordered memory allocator if available.
* If not, use The_Arena().
*/
class PArena
    :
    public Arena
{
public:
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

#ifdef AMREX_CUDA_GE_11_2
private:
    cudaMemPool_t m_pool;
    cuuint64_t m_old_release_threshold;
#endif
};

}
#endif
