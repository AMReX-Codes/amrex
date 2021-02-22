#ifndef AMREX_PARENA_H_
#define AMREX_PARENA_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>

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

#if (__CUDACC_VER_MAJOR__ > 11 || ((__CUDACC_VER_MAJOR__ == 11) && (__CUDACC_VER_MINOR__ >= 2)))
private:
    cudaMemPool_t m_pool;
    cuuint64_t m_old_release_threshold;
#endif
};

}
#endif
