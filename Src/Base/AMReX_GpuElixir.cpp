
#include <AMReX_GpuElixir.H>
#include <AMReX_GpuDevice.H>
#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <memory>

namespace amrex::Gpu {

namespace {

#if defined(AMREX_USE_GPU) && !defined(AMREX_USE_SYCL)

extern "C" {
#if defined(AMREX_USE_HIP)
    void amrex_elixir_delete ( hipStream_t /*stream*/,  hipError_t /*error*/, void* p)
#elif defined(AMREX_USE_CUDA)
    void CUDART_CB amrex_elixir_delete (void* p)
#endif
    {
        auto p_pa = reinterpret_cast<Vector<std::pair<void*,Arena*> >*>(p);
        for (auto const& pa : *p_pa) {
            pa.second->free(pa.first);
        }
        delete p_pa;
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
        if (!m_pa.empty()) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
            auto p = new Vector<std::pair<void*,Arena*> >(std::move(m_pa));
#if defined(AMREX_USE_HIP)
            AMREX_HIP_SAFE_CALL ( hipStreamAddCallback(Gpu::gpuStream(),
                                                       amrex_elixir_delete, (void*)p, 0));
#elif defined(AMREX_USE_CUDA)
            AMREX_CUDA_SAFE_CALL(cudaLaunchHostFunc(Gpu::gpuStream(),
                                                    amrex_elixir_delete, (void*)p));
#endif
#elif defined(AMREX_USE_SYCL)
#ifdef AMREX_USE_CODEPLAY_HOST_TASK
            auto lpa = std::move(m_pa);
            auto& q = *(Gpu::gpuStream().queue);
            try {
                q.submit([&] (sycl::handler& h) {
                    h.codeplay_host_task([=] () {
                        for (auto const& pa : lpa) {
                            pa.second->free(pa.first);
                        }
                    });
                });
            } catch (sycl::exception const& ex) {
                amrex::Abort(std::string("host_task: ")+ex.what()+"!!!!!");
            }
#else
            // xxxxx SYCL todo
            Gpu::streamSynchronize();
            for (auto const& pa : m_pa) {
                pa.second->free(pa.first);
            }
#endif
#endif
        }
    }
    else
#endif
    {
        for (auto const& pa : m_pa) {
            pa.second->free(pa.first);
        }
    }
    m_pa.clear();
}

}
