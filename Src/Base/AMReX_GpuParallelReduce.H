#ifndef AMREX_GPU_PARALLEL_REDUCE_H_
#define AMREX_GPU_PARALLEL_REDUCE_H_
#include <AMReX_Config.H>

#include <AMReX_GpuContainers.H>
#include <AMReX_INT.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <cstddef>

//
// GPU-aware MPI collectives that operate in place on a Gpu::DeviceVector.
//
// These overloads complement the pointer/scalar overloads in
// AMReX_ParallelReduce.H and AMReX_ParallelDescriptor.H. They live in a
// separate header (rather than in those low-level headers) because they need
// the Gpu container/copy machinery (AMReX_GpuContainers.H) that would bloat the
// ParallelReduce.H headers. This header is also pulled into the AMReX_Gpu.H
// umbrella for convenience.
//
// When AMReX is configured with GPU-aware MPI (ParallelDescriptor::UseGpuAwareMpi())
// the device buffer is handed to MPI directly, otherwise the data is staged
// through host (pinned) memory for the collective.
//

namespace amrex {

namespace ParallelAllReduce {

/// \ingroup amrex_mpi
/** Sum-reduce the elements of a Gpu::DeviceVector in place, on all ranks. */
template <typename T>
void Sum (Gpu::DeviceVector<T>& v, MPI_Comm comm)
{
    // GPU-unaware case
#if defined(AMREX_USE_MPI) && defined(AMREX_USE_GPU)
    if (!ParallelDescriptor::UseGpuAwareMpi()) {
        Gpu::PinnedVector<T> hv(v.size());
        Gpu::copy(Gpu::deviceToHost, v.begin(), v.end(), hv.begin());
        Sum(hv.data(), static_cast<int>(hv.size()), comm);
        Gpu::copy(Gpu::hostToDevice, hv.begin(), hv.end(), v.begin());
        return;
    }
#endif

    // GPU-aware case
    Sum(v.data(), static_cast<int>(v.size()), comm);
}

} // namespace ParallelAllReduce

namespace ParallelReduce {

/// \ingroup amrex_mpi
/** Sum-reduce the elements of a Gpu::DeviceVector in place, onto \p root. */
template <typename T>
void Sum (Gpu::DeviceVector<T>& v, int root, MPI_Comm comm)
{
    // GPU-unaware case
#if defined(AMREX_USE_MPI) && defined(AMREX_USE_GPU)
    if (!ParallelDescriptor::UseGpuAwareMpi()) {
        Gpu::PinnedVector<T> hv(v.size());
        // every rank stages its contribution to host for the reduction
        Gpu::copy(Gpu::deviceToHost, v.begin(), v.end(), hv.begin());
        Sum(hv.data(), static_cast<int>(hv.size()), root, comm);
        // only the root receives the reduced result, so only it copies back
        if (ParallelDescriptor::MyProc(comm) == root) {
            Gpu::copy(Gpu::hostToDevice, hv.begin(), hv.end(), v.begin());
        }
        return;
    }
#endif

    // GPU-aware case
    Sum(v.data(), static_cast<int>(v.size()), root, comm);
}

} // namespace ParallelReduce

namespace ParallelDescriptor {

/// \ingroup amrex_mpi
/** Broadcast a Gpu::DeviceVector from \p root to all ranks in \p comm.
 *
 * Contract: on every receiving rank \p v must already be allocated to the same
 * length as on \p root (the caller is responsible for sizing it, e.g. via
 * resize). Calling this with a wrong-sized or empty receiver is undefined
 * behavior. This mirrors the pointer overload Bcast(T*, n, root), which also
 * requires a correctly-sized buffer on every rank.
 *
 * In a debug build (AMREX_DEBUG) the contract is checked: the root's length is
 * broadcast first and asserted against each rank's size before the data are
 * broadcast.
 *
 * The device buffer is passed to MPI directly when GPU-aware MPI is enabled
 * (UseGpuAwareMpi()), otherwise it is staged through host (pinned) memory.
 */
template <typename T>
void Bcast (Gpu::DeviceVector<T>& v, int root, MPI_Comm comm)
{
#ifdef AMREX_USE_MPI
    auto const n = v.size();

#ifdef AMREX_DEBUG
    // verify the pre-allocation contract (the length broadcast happens on every
    // rank, so it is collectively safe and cannot deadlock)
    Long n_root = static_cast<Long>(n);
    Bcast(&n_root, std::size_t(1), root, comm);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_root == static_cast<Long>(n),
        "ParallelDescriptor::Bcast(Gpu::DeviceVector): receiver must be pre-allocated to the root's length");
#endif

    // trivial case: 1 rank
    if (n == 0) { return; }

    // GPU-unaware case
#ifdef AMREX_USE_GPU
    if (!UseGpuAwareMpi()) {
        Gpu::PinnedVector<T> hv(n);
        const bool is_root = (MyProc(comm) == root);
        // only the root needs to stage its data to host before the broadcast
        if (is_root) {
            Gpu::copy(Gpu::deviceToHost, v.begin(), v.end(), hv.begin());
        }
        Bcast(hv.data(), static_cast<std::size_t>(n), root, comm);
        // only the receivers need to copy the broadcast result back to device
        if (!is_root) {
            Gpu::copy(Gpu::hostToDevice, hv.begin(), hv.end(), v.begin());
        }
        return;
    }
#endif

    // GPU-aware case
    Bcast(v.data(), static_cast<std::size_t>(n), root, comm);

#else  // AMREX_USE_MPI
    amrex::ignore_unused(v, root, comm);
#endif
}

} // namespace ParallelDescriptor

} // namespace amrex

#endif /*AMREX_GPU_PARALLEL_REDUCE_H_*/
