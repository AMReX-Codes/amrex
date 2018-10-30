
#ifdef AMREX_USE_CUDA

#include <AMReX_CudaFabImpl.H>
#include <AMReX_CudaDevice.H>
#include <AMReX_Vector.H>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {
namespace Cuda {

namespace {
    bool initialized = false;
    Vector<Vector<std::unique_ptr<FArrayBox> > > fab_stacks;

    inline Vector<std::unique_ptr<FArrayBox> >& get_stack ()
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        return fab_stacks[tid];
    }
}

void
DeviceFabImpl::Initialize ()
{
    if (initialized) return;
    initialized = true;

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif
    fab_stacks.resize(nthreads);
}

void
DeviceFabImpl::Finalize ()
{
    fab_stacks.clear();
}


DeviceFabImpl::DeviceFabImpl ()
{
    static_assert(AMREX_IS_TRIVIALLY_COPYABLE(BaseFabData<Real>),
                  "BaseFabData must be trivially copyable");
    auto& fabstack = get_stack();
    if (fabstack.empty()) {
        m_gpu_fab.reset(new FArrayBox());
    } else {
        m_gpu_fab = std::move(fabstack.back());
        fabstack.pop_back();
        copy_htod();
    }
}

DeviceFabImpl::DeviceFabImpl (Box const& bx, int ncomp)
    : m_cpu_fab(bx,ncomp)
{
    auto& fabstack = get_stack();
    if (fabstack.empty()) {
        m_gpu_fab.reset(new FArrayBox());
    } else {
        m_gpu_fab = std::move(fabstack.back());
        fabstack.pop_back();
    }   
    copy_htod();
}

DeviceFabImpl::~DeviceFabImpl ()
{
    auto& fabstack = get_stack();
    fabstack.push_back(std::move(m_gpu_fab));
}

void
DeviceFabImpl::resize (Box const& bx, int ncomp)
{
    if (ncomp != m_cpu_fab.nComp() || bx != m_cpu_fab.box())
    {
        m_cpu_fab.resize(bx,ncomp);
        copy_htod();
    }
}

FArrayBox*
DeviceFabImpl::fabPtr ()
{
    AMREX_ASSERT(m_gpu_fab != nullptr);
    return m_gpu_fab.get();
}

void
DeviceFabImpl::copy_htod ()
{
    auto dest = static_cast<BaseFabData<Real>*>(m_gpu_fab.get());
    if (inLaunchRegion())
    {
        m_cpu_fab_data = m_cpu_fab;
        m_cpu_fab_data.setOwner(false);
        CudaAPICheck(cudaMemcpyAsync(dest, &m_cpu_fab_data, sizeof(BaseFabData<Real>),
                                     cudaMemcpyHostToDevice, Device::cudaStream()));
    }
    else
    {
        auto src  = static_cast<BaseFabData<Real>*>(&m_cpu_fab);
        std::memcpy(dest, src, sizeof(BaseFabData<Real>));
    }
}

}
}

#else

namespace amrex {
namespace Cuda {

void DeviceFabImpl::Initialize () {}
void DeviceFabImpl::Finalize () {}

DeviceFabImpl::DeviceFabImpl () {}

DeviceFabImpl::DeviceFabImpl (Box const& bx, int ncomp) : m_cpu_fab(bx,ncomp) {}

DeviceFabImpl::~DeviceFabImpl () {}

void DeviceFabImpl::resize (Box const& bx, int ncomp) { m_cpu_fab.resize(bx,ncomp); }

FArrayBox* DeviceFabImpl::fabPtr () { return &m_cpu_fab; }

}}

#endif
