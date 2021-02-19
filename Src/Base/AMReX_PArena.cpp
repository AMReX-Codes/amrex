#include <AMReX_PArena.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuElixir.H>
#include <AMReX_MemPool.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

namespace amrex {

PArena::PArena (Long release_threshold)
{
#if (__CUDACC_VER_MAJOR__ > 11 || ((__CUDACC_VER_MAJOR__ == 11) && (__CUDACC_VER_MINOR__ >= 2)))
    AMREX_CUDA_SAFE_CALL(cudaDeviceGetMemPool(&m_pool, Gpu::Device::deviceId()));
    AMREX_CUDA_SAFE_CALL(cudaMemPoolGetAttribute(m_pool, cudaMemPoolAttrReleaseThreshold,
                                                 &m_old_release_threshold));
    cuuint64_t value = release_threshold;
    AMREX_CUDA_SAFE_CALL(cudaMemPoolSetAttribute(m_pool, cudaMemPoolAttrReleaseThreshold, &value));
#else
    amrex::ignore_unused(release_threshold);
#endif
}

PArena::~PArena ()
{
#if (__CUDACC_VER_MAJOR__ > 11 || ((__CUDACC_VER_MAJOR__ == 11) && (__CUDACC_VER_MINOR__ >= 2)))
    AMREX_CUDA_SAFE_CALL(cudaMemPoolSetAttribute(m_pool, cudaMemPoolAttrReleaseThreshold,
                                                 &m_old_release_threshold));
#endif
}

void*
PArena::alloc (std::size_t nbytes)
{
#if defined(AMREX_USE_GPU)

#if (__CUDACC_VER_MAJOR__ > 11 || ((__CUDACC_VER_MAJOR__ == 11) && (__CUDACC_VER_MINOR__ >= 2)))
    if (Gpu::Device::memoryPoolsSupported()) {
        void* p;
        AMREX_CUDA_SAFE_CALL(cudaMallocAsync(&p, nbytes, m_pool, Gpu::gpuStream()));
        return p;
    } else
#endif
    {
        return The_Arena()->alloc(nbytes);
    }

#elif defined(AMREX_USE_OMP)

    if (omp_in_parallel()) {
        return amrex_mempool_alloc(nbytes);
    } else {
        return The_Arena()->alloc(nbytes);
    }

#else

    return The_Arena()->alloc(nbytes);

#endif
}

void
PArena::free (void* p)
{
#if defined(AMREX_USE_GPU)

#if (__CUDACC_VER_MAJOR__ > 11 || ((__CUDACC_VER_MAJOR__ == 11) && (__CUDACC_VER_MINOR__ >= 2)) )
    if (Gpu::Device::memoryPoolsSupported()) {
        AMREX_CUDA_SAFE_CALL(cudaFreeAsync(p, Gpu::gpuStream()));
    } else
#endif
    {
        Elixir eli(p, The_Arena());
    }

#elif defined(AMREX_USE_OMP)

    if (omp_in_parallel()) {
        amrex_mempool_free(p);
    } else {
        The_Arena()->free(p);
    }

#else

    The_Arena()->free(p);

#endif
}

}
