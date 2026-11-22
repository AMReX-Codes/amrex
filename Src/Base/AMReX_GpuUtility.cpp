
#include <AMReX_GpuUtility.H>

#include <atomic>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

namespace amrex { // NOLINT(modernize-concat-nested-namespaces)

#ifdef AMREX_USE_GPU

std::ostream&
operator<< (std::ostream& os, const dim3& d)
{
    os << '(' << d.x << ' ' << d.y << ' ' << d.z << ')';
    return os;
}

#endif

namespace Gpu {

#ifdef AMREX_USE_GPU
namespace {
    std::atomic<int> s_stream_iter_depth{0};
}
#endif

int
StreamIter::currentDepth () noexcept
{
#ifdef AMREX_USE_GPU
    return s_stream_iter_depth.load(std::memory_order_relaxed);
#else
    return 0;
#endif
}

StreamIter::StreamIter (int n, bool is_thread_safe) noexcept
    : m_n(n), m_i(0), m_threadsafe(is_thread_safe), m_sync(!Gpu::inStreamOrderedRegion())
{
    init();
}

StreamIter::StreamIter (int n, const StreamItInfo& info, bool is_thread_safe) noexcept
    : m_n(n), m_i(0), m_threadsafe(is_thread_safe), m_sync(info.device_sync)
{
    init();
}

void
StreamIter::init () noexcept // NOLINT
{
    amrex::ignore_unused(m_threadsafe);
    amrex::ignore_unused(m_sync, m_join);
#if defined(AMREX_USE_GPU)
    ++s_stream_iter_depth;
    if (m_sync) {
#ifdef AMREX_USE_OMP
#pragma omp single
#endif
        Gpu::streamSynchronize();
    } else if (Gpu::inStreamOrderedRegion() && Gpu::numGpuStreams() > 1) {
        m_join = true;
        m_origin = Gpu::Device::streamIndex();
#ifdef AMREX_USE_OMP
#pragma omp single
#endif
        Gpu::Device::forkStreams(m_n, m_origin);
    }
    Gpu::Device::setStreamIndex(m_i);
#elif defined(AMREX_USE_OMP)
    int nthreads = omp_get_num_threads();
    if (nthreads > 1) {
        int tid = omp_get_thread_num();
        int nr = m_n / nthreads;
        int nlft = m_n - nr*nthreads;
        if (tid < nlft) { // get nr+1 items
            m_i = tid * (nr+1);
            m_n = m_i + nr+1;
        } else {         // get nr items
            m_i = tid * nr + nlft;
            m_n = m_i + nr;
        }
    }
#endif
}

StreamIter::~StreamIter () { // NOLINT(modernize-use-equals-default)
#ifdef AMREX_USE_GPU
    if (m_sync) {
        const int nstreams = std::min(m_n, Gpu::numGpuStreams());
        for (int i = 0; i < nstreams; ++i) {
            Gpu::Device::setStreamIndex(i);
            Gpu::streamSynchronize();
        }
    } else if (m_join) {
        Gpu::Device::joinStreams(m_n, m_origin);
    }
    AMREX_GPU_ERROR_CHECK();
    if (m_join) {
        Gpu::Device::setStreamIndex(m_origin);
    } else {
        Gpu::Device::resetStreamIndex();
    }
    --s_stream_iter_depth;
#endif
}

#ifdef AMREX_USE_GPU
void
StreamIter::operator++ ()
{
    ++m_i;
    if (m_threadsafe) {
        Gpu::Device::setStreamIndex(m_i);
        AMREX_GPU_ERROR_CHECK();
    }
}
#endif

}}
