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
    PArena& operator= (const PArena& rhs) = delete;
    virtual ~PArena () override;

    virtual void* alloc (std::size_t nbytes) override final;
    virtual void free (void* p) override final;

    virtual bool isDeviceAccessible () const override final;
    virtual bool isHostAccessible () const override final;

    virtual bool isManaged () const override final;
    virtual bool isDevice () const override final;
    virtual bool isPinned () const override final;

#ifdef AMREX_CUDA_GE_11_2
private:
    cudaMemPool_t m_pool;
    cuuint64_t m_old_release_threshold;
#endif
};

}
#endif
