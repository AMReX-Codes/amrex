
#include <AMReX_GpuUtility.H>

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

StreamIter::StreamIter (int n, bool is_thread_safe) noexcept
    // Synchronizes unless in a Gpu::StreamOrderedRegion whose no-sync mode
    // is still on (Gpu::setNoSyncRegion(false) turns host syncs back on).
    : m_n(n), m_i(0), m_threadsafe(is_thread_safe),
      m_sync(!(Gpu::inStreamOrderedRegion() && Gpu::inNoSyncRegion()))
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
    amrex::ignore_unused(m_threadsafe, m_sync);
#if defined(AMREX_USE_GPU)
    // In a Gpu::StreamOrderedRegion, a loop that does not synchronize is
    // ordered on the GPU instead. Outside it, such a loop is not ordered. A
    // loop that is not thread-safe uses one stream.
    const bool ordered = Gpu::inStreamOrderedRegion();
    const auto edge = m_sync ? detail::StreamEdge::Host
                             : (ordered ? detail::StreamEdge::Gpu : detail::StreamEdge::None);
#ifdef AMREX_USE_OMP
    const bool in_team = omp_in_parallel() && omp_get_num_threads() > 1;
#else
    const bool in_team = false;
#endif
    detail::theStreamSequencer().fork(m_loop_streams, {.nitems = m_n,
                                        .max_width = m_threadsafe ? Gpu::numGpuStreams() : 1,
                                        .fork = edge, .join = edge, .shared = in_team,
                                        .ordered = ordered});
    if (m_n > 0) { detail::theStreamSequencer().lane(m_loop_streams, m_i); }
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
    detail::theStreamSequencer().join(m_loop_streams);
    AMREX_GPU_ERROR_CHECK();
#endif
}

#ifdef AMREX_USE_GPU
void
StreamIter::operator++ ()
{
    ++m_i;
    if (m_i < m_n) { detail::theStreamSequencer().lane(m_loop_streams, m_i); }
    if (m_threadsafe) { AMREX_GPU_ERROR_CHECK(); }
}
#endif

}}
