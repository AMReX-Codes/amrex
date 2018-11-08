
#include <AMReX.H>
#include <AMReX_CudaAsyncFab.H>
#include <AMReX_CudaAsyncFabImpl.H>
#include <AMReX_CudaDevice.H>
#include <mutex>

#ifdef AMREX_USE_CUDA

namespace {
    std::mutex cudafab_callback_mutex;
}

extern "C" {
    void CUDART_CB amrex_devicefab_delete (cudaStream_t stream, cudaError_t error, void* p) {
        std::lock_guard<std::mutex> guard(cudafab_callback_mutex);
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
{
    m_impl.reset(new AsyncFabImpl());
}

AsyncFab::AsyncFab (Box const& bx, int ncomp)
{
    m_impl.reset(new AsyncFabImpl(bx,ncomp));
}

AsyncFab::AsyncFab (FArrayBox& a_fab)
{
    m_impl.reset(new AsyncFabImpl(a_fab));
}

AsyncFab::AsyncFab (FArrayBox& a_fab, Box const& bx, int ncomp)
{
    m_impl.reset(new AsyncFabImpl(a_fab, bx, ncomp));
}

void
AsyncFab::clear ()
{
#ifdef AMREX_USE_CUDA
    if (inLaunchRegion())
    {
        if (m_impl != nullptr) {
            AsyncFabImpl* p = m_impl.release();
// CUDA 10        CudaAPICheck(cudaLaunchHostFunc(Device::cudaStream(), amrex_devicefab_delete, p));
            CudaAPICheck(cudaStreamAddCallback(Device::cudaStream(), amrex_devicefab_delete, p, 0));
        }
    }
    else
#endif
    {
        m_impl.reset();
    }
}

FArrayBox*
AsyncFab::fabPtr ()
{
    return m_impl->fabPtr();
}

AsyncFab::~AsyncFab ()
{
    clear();
}

}
}
