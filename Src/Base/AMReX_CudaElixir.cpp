
#include <AMReX_CudaElixir.H>
#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <AMReX_CudaDevice.H>

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
    if (inLaunchRegion())
    {
        if (m_p != nullptr) {
            void** p = static_cast<void**>(std::malloc(2*sizeof(void*)));
            p[0] = m_p;
            p[1] = (void*)m_arena;
            AMREX_GPU_SAFE_CALL(cudaStreamAddCallback(Device::cudaStream(),
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
