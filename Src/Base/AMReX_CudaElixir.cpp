
#include <AMReX_CudaElixir.H>
#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <AMReX_GpuDevice.H>

namespace amrex {
namespace Cuda {

namespace {

#ifdef AMREX_USE_CUDA

extern "C" {
    void CUDART_CB amrex_elixir_delete (cudaStream_t stream, cudaError_t error, void* p)
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
#ifdef AMREX_USE_CUDA
    if (Gpu::inLaunchRegion())
    {
        if (m_p != nullptr) {
            void** p = static_cast<void**>(std::malloc(2*sizeof(void*)));
            p[0] = m_p;
            p[1] = (void*)m_arena;
            AMREX_CUDA_SAFE_CALL(cudaStreamAddCallback(Gpu::gpuStream(),
                                                       amrex_elixir_delete, p, 0));
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
