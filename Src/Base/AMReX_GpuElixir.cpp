
#include <AMReX_GpuElixir.H>
#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <AMReX_GpuDevice.H>

namespace amrex {
namespace Gpu {

namespace {

#if defined(AMREX_USE_GPU) && !defined(AMREX_USE_DPCPP)

extern "C" {
#if defined(AMREX_USE_HIP)
    void amrex_elixir_delete ( hipStream_t /*stream*/,  hipError_t /*error*/, void* p)
#elif defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10)
    void CUDART_CB amrex_elixir_delete (void* p)
#elif defined(AMREX_USE_CUDA)
    void CUDART_CB amrex_elixir_delete (cudaStream_t /*stream*/, cudaError_t /*error*/, void* p)
#endif
    {
        void** pp = (void**)p;
        void* d = pp[0];
        Arena* arena = (Arena*)(pp[1]);
        std::free(p);
        arena->free(d);
    }
}

#endif

}

void
Elixir::clear () noexcept
{
#if defined(AMREX_USE_GPU)
    if (Gpu::inLaunchRegion())
    {
        if (m_p != nullptr) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
            void** p = static_cast<void**>(std::malloc(2*sizeof(void*)));
            p[0] = m_p;
            p[1] = (void*)m_arena;
#if defined(AMREX_USE_HIP)
            AMREX_HIP_SAFE_CALL ( hipStreamAddCallback(Gpu::gpuStream(),
                                                       amrex_elixir_delete, p, 0));
#elif defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10)
            AMREX_CUDA_SAFE_CALL(cudaLaunchHostFunc(Gpu::gpuStream(),
                                                    amrex_elixir_delete, p));
#elif defined(AMREX_USE_CUDA)
            AMREX_CUDA_SAFE_CALL(cudaStreamAddCallback(Gpu::gpuStream(),
                                                       amrex_elixir_delete, p, 0));
#endif
            Gpu::callbackAdded();
#elif defined(AMREX_USE_DPCPP)
            // xxxxx DPCPP todo
            Gpu::streamSynchronize();
            if (m_p != nullptr) m_arena->free(m_p);
#endif
        }
    }
    else
#endif
    {
        if (m_p != nullptr) m_arena->free(m_p);
    }
    m_p = nullptr;
    m_arena = nullptr;
}

}
}
