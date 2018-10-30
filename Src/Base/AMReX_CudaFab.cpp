
#include <AMReX.H>
#include <AMReX_CudaFab.H>
#include <AMReX_CudaFabImpl.H>
#include <AMReX_CudaDevice.H>
#include <AMReX_TinyProfiler.H>

#ifdef AMREX_USE_CUDA
extern "C" {
    void CUDART_CB amrex_devicefab_delete (cudaStream_t stream, cudaError_t error, void* p) {
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

void
DeviceFab::resize (Box const& bx, int ncomp)
{
    m_impl->resize(bx, ncomp);
}

FArrayBox*
DeviceFab::fabPtr ()
{
    return m_impl->fabPtr();
}

DeviceFab::~DeviceFab ()
{
#ifdef AMREX_USE_CUDA
    BL_PROFILE("DeviceFab::~DeviceFab()");
    if (inLaunchRegion())
    {
        DeviceFabImpl* p = m_impl.release();
// CUDA 10        CudaAPICheck(cudaLaunchHostFunc(Device::cudaStream(), amrex_devicefab_delete, p));
        CudaAPICheck(cudaStreamAddCallback(Device::cudaStream(), amrex_devicefab_delete, p, 0));
    }
    else
#endif
    {
        m_impl.reset();
    }
}

}
}
