
#include <AMReX.H>
#include <AMReX_CudaFab.H>
#include <AMReX_CudaFabImpl.H>
#include <AMReX_CudaDevice.H>
#include <mutex>

#ifdef AMREX_USE_CUDA

namespace {
    std::mutex callback_mutex;
}

extern "C" {
    void CUDART_CB amrex_devicefab_delete (cudaStream_t stream, cudaError_t error, void* p) {
        std::lock_guard<std::mutex> guard(callback_mutex);
        delete (amrex::Cuda::DeviceFabImpl*)p;
    }
}

#endif

namespace amrex {
namespace Cuda {

void
DeviceFab::Initialize ()
{
    DeviceFabImpl::Initialize();
    amrex::ExecOnFinalize(DeviceFab::Finalize);
}

void
DeviceFab::Finalize ()
{
    DeviceFabImpl::Finalize();
}

DeviceFab::DeviceFab ()
{
    m_impl.reset(new DeviceFabImpl());
}

DeviceFab::DeviceFab (Box const& bx, int ncomp)
{
    m_impl.reset(new DeviceFabImpl(bx,ncomp));
}

DeviceFab::DeviceFab (FArrayBox& a_fab)
{
    m_impl.reset(new DeviceFabImpl(a_fab));
}

DeviceFab::DeviceFab (FArrayBox& a_fab, Box const& bx, int ncomp)
{
    m_impl.reset(new DeviceFabImpl(a_fab, bx, ncomp));
}

void
DeviceFab::clear ()
{
#ifdef AMREX_USE_CUDA
    if (inLaunchRegion())
    {
        if (m_impl != nullptr) {
            DeviceFabImpl* p = m_impl.release();
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
DeviceFab::fabPtr ()
{
    return m_impl->fabPtr();
}

DeviceFab::~DeviceFab ()
{
    clear();
}

}
}
