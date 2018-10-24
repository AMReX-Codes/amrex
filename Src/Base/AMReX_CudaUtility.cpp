
#include <AMReX_CudaUtility.H>
#include <AMReX_Device.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

#if AMREX_USE_CUDA

std::ostream&
operator<< (std::ostream& os, const dim3& d)
{
    os << '(' << d.x << ' ' << d.y << ' ' << d.z << ')';
    return os;
}

#endif

namespace Cuda {

StreamIter::StreamIter (const int n)
    : m_n(n), m_i(0)
{
#if defined(AMREX_USE_CUDA)
    Device::set_stream_index(m_i);
#elif defined(_OPENMP)
    int nthreads = omp_get_num_threads();
    if (nthreads > 1) {
        int tid = omp_get_thread_num();
        int nr = n / nthreads;
        int nlft = n - nr*nthreads;
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

StreamIter::~StreamIter () {
#ifdef AMREX_USE_CUDA
    Device::synchronize();
    Device::check_for_errors();
    Device::set_stream_index(-1);
#endif
}

#ifdef AMREX_USE_CUDA
void
StreamIter::operator++ ()
{
    ++m_i;
    Device::set_stream_index(m_i);
    Device::check_for_errors();
}
#endif

}}

