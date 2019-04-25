
#include <AMReX.H>
#include <AMReX_CudaAsyncFab.H>
#include <AMReX_CudaAsyncFabImpl.H>
#include <AMReX_CudaDevice.H>

#ifdef AMREX_USE_CUDA

extern "C" {
    void CUDART_CB amrex_devicefab_delete (cudaStream_t stream, cudaError_t error, void* p) {
        delete (amrex::Cuda::AsyncFabImpl*)p;
    }
}

#endif

namespace amrex {
namespace Cuda {

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
#ifdef AMREX_USE_CUDA
    if (inLaunchRegion())
    {
        if (m_impl != nullptr) {
// CUDA 10        AMREX_GPU_SAFE_CALL(cudaLaunchHostFunc(Device::cudaStream(), amrex_devicefab_delete, p));
            AMREX_GPU_SAFE_CALL(cudaStreamAddCallback(Device::cudaStream(),
                                                      amrex_devicefab_delete,
                                                      m_impl, 0));
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
