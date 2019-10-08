
#include <AMReX_GpuElixir.H>
#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <AMReX_GpuDevice.H>

namespace amrex {
namespace Gpu {

namespace {

#ifdef AMREX_USE_GPU

extern "C" {
AMREX_HIP_OR_CUDA(
    void  HIPRT_CB amrex_elixir_delete ( hipStream_t stream,  hipError_t error, void* p),
    void CUDART_CB amrex_elixir_delete (cudaStream_t stream, cudaError_t error, void* p))
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
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        if (m_p != nullptr) {
            void** p = static_cast<void**>(std::malloc(2*sizeof(void*)));
            p[0] = m_p;
            p[1] = (void*)m_arena;
            AMREX_HIP_OR_CUDA(
                AMREX_HIP_SAFE_CALL ( hipStreamAddCallback(Gpu::gpuStream(),
                                                           amrex_elixir_delete, p, 0));,
                AMREX_CUDA_SAFE_CALL(cudaStreamAddCallback(Gpu::gpuStream(),
                                                           amrex_elixir_delete, p, 0)););
            Gpu::callbackAdded();
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
