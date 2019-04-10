#include <AMReX_CudaAsyncFabImpl.H>

#ifdef AMREX_USE_CUDA

#include <AMReX_CudaDevice.H>
#include <AMReX_Vector.H>
#include <cstring>
#include <mutex>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {
namespace Cuda {

#ifdef AMREX_FAB_IS_PINNED

void AsyncFabImpl::Initialize () {}
void AsyncFabImpl::Finalize () {}

AsyncFabImpl::AsyncFabImpl () : m_gpu_fab(new FArrayBox()) {}

AsyncFabImpl::AsyncFabImpl (Box const& bx, int ncomp) : m_gpu_fab(new FArrayBox(bx,ncomp)) {}

AsyncFabImpl::AsyncFabImpl (FArrayBox& a_fab)
    : m_gpu_fab(a_fab.isAllocated()
                ? new FArrayBox(a_fab.box(), a_fab.nComp())
                : new FArrayBox())
{}

AsyncFabImpl::AsyncFabImpl (FArrayBox& /*a_fab*/, Box const& bx, int ncomp)
    : AsyncFabImpl(bx,ncomp)
{}

AsyncFabImpl::~AsyncFabImpl () {}

FArrayBox*
AsyncFabImpl::fabPtr () noexcept
{
    AMREX_ASSERT(m_gpu_fab != nullptr);
    return m_gpu_fab.get();
}

FArrayBox& AsyncFabImpl::hostFab () noexcept { return *m_gpu_fab; }

#else

//
// Fab is Managed
//

namespace {
    bool initialized = false;

    std::mutex asyncfab_mutex;

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
AsyncFabImpl::Initialize ()
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
AsyncFabImpl::Finalize ()
{
    fab_stacks.clear();
}


AsyncFabImpl::AsyncFabImpl ()
{
    static_assert(AMREX_IS_TRIVIALLY_COPYABLE(BaseFabData<Real>),
                  "BaseFabData must be trivially copyable");
    auto& fabstack = get_stack();
    bool do_copy = false;
    {
        std::lock_guard<std::mutex> guard(asyncfab_mutex);
        if (fabstack.empty()) {
            m_gpu_fab.reset(new FArrayBox());
        } else {
            m_gpu_fab = std::move(fabstack.back());
            fabstack.pop_back();
            do_copy = true;
        }
    }
    if (do_copy) copy_htod();
}

AsyncFabImpl::AsyncFabImpl (Box const& bx, int ncomp)
    : m_cpu_fab(bx,ncomp)
{
    auto& fabstack = get_stack();
    {
        std::lock_guard<std::mutex> guard(asyncfab_mutex);
        if (fabstack.empty()) {
            m_gpu_fab.reset(new FArrayBox());  // yes, we build an empty fab here, later it will be overwritten by copy_htod
        } else {
            m_gpu_fab = std::move(fabstack.back());
            fabstack.pop_back();
        }
    }
    copy_htod();
}

AsyncFabImpl::AsyncFabImpl (FArrayBox& a_fab)
{
    if (a_fab.isAllocated()) {
        m_cpu_fab.resize(a_fab.box(), a_fab.nComp());
    }
    auto& fabstack = get_stack();
    {
        std::lock_guard<std::mutex> guard(asyncfab_mutex);
        if (fabstack.empty()) {
            m_gpu_fab.reset(new FArrayBox());  // yes, we build an empty fab here, later it will be overwritten by copy_htod
        } else {
            m_gpu_fab = std::move(fabstack.back());
            fabstack.pop_back();
        }
    }
    copy_htod();
}

AsyncFabImpl::AsyncFabImpl (FArrayBox& /*a_fab*/, Box const& bx, int ncomp)
    : AsyncFabImpl(bx,ncomp)
{}

AsyncFabImpl::~AsyncFabImpl ()
{
    std::lock_guard<std::mutex> guard(asyncfab_mutex);
    auto& fabstack = get_stack();
    fabstack.push_back(std::move(m_gpu_fab));
}

FArrayBox*
AsyncFabImpl::fabPtr () noexcept
{
    AMREX_ASSERT(m_gpu_fab != nullptr);
    return m_gpu_fab.get();
}

FArrayBox&
AsyncFabImpl::hostFab () noexcept
{
    return m_cpu_fab;
}

void
AsyncFabImpl::copy_htod ()
{
    auto dest = static_cast<BaseFabData<Real>*>(m_gpu_fab.get());
    if (inLaunchRegion())
    {
        m_cpu_fab_data = m_cpu_fab;
        m_cpu_fab_data.setOwner(false);
        AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(dest, &m_cpu_fab_data, sizeof(BaseFabData<Real>),
                                            cudaMemcpyHostToDevice, Device::cudaStream()));
    }
    else
    {
        auto src  = static_cast<BaseFabData<Real>*>(&m_cpu_fab);
        std::memcpy(dest, src, sizeof(BaseFabData<Real>));
        dest->setOwner(false);
    }
}

#endif

}
}

#else

namespace amrex {
namespace Cuda {

void AsyncFabImpl::Initialize () {}
void AsyncFabImpl::Finalize () {}

AsyncFabImpl::AsyncFabImpl () {}

AsyncFabImpl::AsyncFabImpl (Box const& bx, int ncomp)
    : m_cpu_fab(bx,ncomp), m_cpu_fab_alias(&m_cpu_fab) {}

AsyncFabImpl::AsyncFabImpl (FArrayBox& a_fab) : m_cpu_fab_alias (&a_fab) {}

AsyncFabImpl::AsyncFabImpl (FArrayBox& a_fab, Box const& bx, int ncomp)
    : m_cpu_fab_alias(&a_fab)
{
    m_cpu_fab_alias->resize(bx,ncomp);
}

AsyncFabImpl::~AsyncFabImpl () {}

FArrayBox* AsyncFabImpl::fabPtr () noexcept { return m_cpu_fab_alias; }

FArrayBox& AsyncFabImpl::hostFab () noexcept { return *m_cpu_fab_alias; }

}}

#endif
