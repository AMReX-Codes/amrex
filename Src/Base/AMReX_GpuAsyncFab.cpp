
#include <AMReX.H>
#include <AMReX_GpuAsyncFab.H>
#include <AMReX_GpuAsyncFabImpl.H>
#include <AMReX_GpuDevice.H>

#ifdef AMREX_USE_GPU

extern "C" {

#if defined(AMREX_USE_HIP) || (__CUDACC_VER_MAJOR__ < 10)
AMREX_HIP_OR_CUDA(
    void  HIPRT_CB amrex_devicefab_delete ( hipStream_t stream,  hipError_t error, void* p),
    void CUDART_CB amrex_devicefab_delete (cudaStream_t stream, cudaError_t error, void* p))
    {
        delete (amrex::Gpu::AsyncFabImpl*)p;
    }
#else
    void CUDART_CB amrex_devicefab_delete (void* p) {
        delete (amrex::Gpu::AsyncFabImpl*)p;
    }
#endif

}

#endif

namespace amrex {
namespace Gpu {

void
AsyncFab::Initialize ()
{
    AsyncFabImpl::Initialize();
    amrex::ExecOnFinalize(AsyncFab::Finalize);
}

void
AsyncFab::Finalize ()
{
    AsyncFabImpl::Finalize();
}

AsyncFab::AsyncFab ()
    : m_impl(new AsyncFabImpl()),
      m_fab(m_impl->fabPtr())
{}

AsyncFab::AsyncFab (Box const& bx, int ncomp)
    : m_impl(new AsyncFabImpl(bx,ncomp)),
      m_fab(m_impl->fabPtr())
{}

AsyncFab::AsyncFab (FArrayBox& a_fab)
    : m_impl(new AsyncFabImpl(a_fab)),
      m_fab(m_impl->fabPtr())
{}

AsyncFab::AsyncFab (FArrayBox& a_fab, Box const& bx, int ncomp)
    : m_impl(new AsyncFabImpl(a_fab, bx, ncomp)),
      m_fab(m_impl->fabPtr())
{}

void
AsyncFab::clear ()
{
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        if (m_impl != nullptr) {
#if defined(AMREX_USE_HIP) || (__CUDACC_VER_MAJOR__ < 10)
            AMREX_HIP_OR_CUDA(
                AMREX_HIP_SAFE_CALL(  hipStreamAddCallback(Gpu::gpuStream(),
                                                           amrex_devicefab_delete,
                                                           m_impl, 0));,
                AMREX_CUDA_SAFE_CALL(cudaStreamAddCallback(Gpu::gpuStream(),
                                                           amrex_devicefab_delete,
                                                           m_impl, 0)););
#else
            cudaHostFn_t clear = amrex_devicefab_delete;
            AMREX_GPU_SAFE_CALL(cudaLaunchHostFunc(Gpu::gpuStream(), clear, m_impl));
#endif
        }
    }
    else
#endif
    {
        delete m_impl;
    }
    m_impl = nullptr;
    m_fab = nullptr;
}

FArrayBox&
AsyncFab::hostFab () const noexcept
{
    return m_impl->hostFab();
}

Array4<Real const>
AsyncFab::array () const noexcept
{
    return m_impl->hostFab().const_array();
}

Array4<Real>
AsyncFab::array () noexcept
{
    return m_impl->hostFab().array();
}

}
}
