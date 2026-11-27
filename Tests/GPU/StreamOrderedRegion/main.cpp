#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_AsyncOut.H>
#include <AMReX_FileSystem.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuUtility.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_OpenMP.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <atomic>
#include <chrono>
#include <map>
#include <mutex>
#include <vector>

using namespace amrex;

/**
 * Tests of the stream ordering of multi-stream loops in a
 * Gpu::StreamOrderedRegion.
 *
 * 1. A loop started while a stream other than stream 0 is current (here,
 *    nested in an MFIter iteration on stream 1) must see the work queued on
 *    that stream before it, and the work queued on that stream after it must
 *    see the loop's work.
 *
 *    Also for loops on a single stream (MFItInfo::SetNumStreams(1), or a
 *    StreamIter in a Gpu::SingleStreamRegion).
 *
 *    A nested loop that does not fork onto other streams (one outside the
 *    launch region, an AllBoxes MFIter, or a StreamIter with device_sync)
 *    must leave the current stream as it found it, so that the work queued
 *    after it on that stream sees the work queued before it.
 *
 * 2. FabArray::arrays() built for the first time in a StreamIter or MFIter
 *    iteration must be usable from the loop's other streams, also after
 *    another loop nested in that iteration has finished.
 *
 * 3. VisMF::AsyncWrite (amrex.async_out=1) in the region writes the data
 *    that the GPU work queued before it produces (smoke test).
 *
 * 4. (With OpenMP.) Memory released with Gpu::freeAsync on a stream while
 *    another thread joins that stream (Gpu::StreamSequencer::join) must not
 *    be freed before the kernels queued ahead of it on that stream finish.
 *    This is a stress test: the window it probes is narrow, and on CUDA it
 *    is usually shorter than a kernel launch, so it mainly catches
 *    regressions on HIP.
 *
 * 5. With OpenMP, an MFIter nested in each iteration of an MFIter run by an
 *    active OpenMP team outside the launch region must not hang when the
 *    team's threads run different numbers of outer iterations.
 *
 * 6. With CUDA, a loop over fewer boxes than there are GPU streams that
 *    builds FabArray::arrays() while the current stream is captured into a
 *    CUDA graph must not pull the streams it does not use into the capture,
 *    and the graph must read the arrays correctly. Also after a loop that
 *    synchronizes (MFItInfo::SetDeviceSync(true)) in the same region.
 *
 * 7. FabArray::arrays() built in a Gpu::StreamOrderedRegion nested in an
 *    iteration of a synchronizing multi-stream MFIter must be usable from
 *    that loop's other streams.
 *
 * 8. An MFIter with no local items keeps the current stream.
 */

namespace {

// Dependent integer arithmetic that keeps one GPU thread busy for a while.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
unsigned spin (Long n) noexcept
{
    unsigned x = 1u;
    for (Long i = 0; i < n; ++i) { x = x * 1664525u + 1013904223u; }
    return x;
}

// Number of spin iterations that take about one millisecond on the GPU.
Long spin_iters_per_ms ()
{
    Gpu::DeviceVector<unsigned> sink(1);
    auto* p = sink.data();
    Long n = 1 << 16;
    while (true) {
        Gpu::synchronize();
        const double t0 = amrex::second();
        amrex::single_task([=] AMREX_GPU_DEVICE () { p[0] = spin(n); });
        Gpu::synchronize();
        const double dt = amrex::second() - t0;
        if (dt > 0.02 || n > (Long(1) << 40)) {
            return std::max(Long(1), static_cast<Long>(static_cast<double>(n) * 1.e-3 / dt));
        }
        n *= 4;
    }
}

enum struct LoopKind {
    MFIter, StreamIterInfo, StreamIterPlain, MFIterOneStream, StreamIterSingleStreamRegion
};
constexpr LoopKind all_loop_kinds[] = {
    LoopKind::MFIter, LoopKind::StreamIterInfo, LoopKind::StreamIterPlain,
    LoopKind::MFIterOneStream, LoopKind::StreamIterSingleStreamRegion
};

const char* loop_name (LoopKind kind)
{
    switch (kind) {
    case LoopKind::MFIter: return "MFIter";
    case LoopKind::StreamIterInfo: return "StreamIter(n, StreamItInfo)";
    case LoopKind::StreamIterPlain: return "StreamIter(n)";
    case LoopKind::MFIterOneStream: return "MFIter(SetNumStreams(1))";
    case LoopKind::StreamIterSingleStreamRegion: return "StreamIter in SingleStreamRegion";
    }
    return "";
}

// Runs a multi-stream loop (MFIter or StreamIter) over the boxes of mf from
// inside an MFIter iteration on stream 1, checks that the loop sees the work
// queued on stream 1 before it, and that work queued after it sees the loop's.
int test_nested_loop (MultiFab const& mf, Long iters_per_ms, LoopKind kind,
                      bool warmup = false)
{
    const char* name = loop_name(kind);
    const int nboxes = mf.local_size();

    Gpu::DeviceVector<int> written(1, 0);  // set on stream 1 before the loop
    Gpu::DeviceVector<int> seen(nboxes, -1);  // what each loop item read
    Gpu::DeviceVector<unsigned> sink(nboxes, 0);
    Gpu::DeviceVector<int> joined(1, -1);  // read on stream 1 after the loop
    auto* pw = written.data();
    auto* ps = seen.data();
    auto* pk = sink.data();
    auto* pj = joined.data();
    const Long long_wait = warmup ? 0 : 50 * iters_per_ms;
    const Long short_wait = warmup ? 0 : iters_per_ms / 4;
    int origin = -1;
    int after = -1;
    Gpu::synchronize();

    {
        Gpu::StreamOrderedRegion stream_ordered;
        for (MFIter omfi(mf); omfi.isValid(); ++omfi) {
            if (omfi.LocalIndex() != 1) { continue; }
            origin = Gpu::Device::streamIndex();

            // Slow producer on the current stream.
            amrex::single_task([=] AMREX_GPU_DEVICE () {
                pk[0] = spin(long_wait);
                pw[0] = 1;
            });

            auto body = [&] (int b) {
                amrex::single_task([=] AMREX_GPU_DEVICE () {
                    const int v = pw[0];
                    pk[b] = spin(short_wait);
                    ps[b] = v;
                });
            };
            if (kind == LoopKind::StreamIterInfo) {
                for (Gpu::StreamIter sit(nboxes, Gpu::StreamItInfo{}); sit.isValid(); ++sit) {
                    body(sit());
                }
            } else if (kind == LoopKind::StreamIterPlain) {
                for (Gpu::StreamIter sit(nboxes); sit.isValid(); ++sit) {
                    body(sit());
                }
            } else if (kind == LoopKind::MFIterOneStream) {
                for (MFIter mfi(mf, MFItInfo().SetNumStreams(1)); mfi.isValid(); ++mfi) {
                    body(mfi.LocalIndex());
                }
            } else if (kind == LoopKind::StreamIterSingleStreamRegion) {
                Gpu::SingleStreamRegion single_stream;
                for (Gpu::StreamIter sit(nboxes); sit.isValid(); ++sit) {
                    body(sit());
                }
            } else {
                for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                    body(mfi.LocalIndex());
                }
            }

            after = Gpu::Device::streamIndex();
            // Consumer of every loop item's result.
            amrex::single_task([=] AMREX_GPU_DEVICE () {
                int s = 0;
                for (int b = 0; b < nboxes; ++b) { s += ps[b]; }
                pj[0] = s;
            });
        }
    }
    Gpu::synchronize();
    if (warmup) { return 0; }

    Vector<int> h_seen(nboxes);
    Gpu::copy(Gpu::deviceToHost, seen.begin(), seen.end(), h_seen.begin());
    int h_joined = 0;
    Gpu::copy(Gpu::deviceToHost, joined.begin(), joined.end(), &h_joined);

    int nfail = 0;
    int nstale = 0;
    for (int v : h_seen) { if (v != 1) { ++nstale; } }
    if (nstale > 0) {
        amrex::Print() << "FAIL nested " << name << ": " << nstale << " of " << nboxes
                       << " loop items ran before the work queued on stream " << origin
                       << " ahead of the loop\n";
        ++nfail;
    }
    if (after != origin) {
        amrex::Print() << "FAIL nested " << name << ": current stream after the loop is "
                       << after << ", expected " << origin << "\n";
        ++nfail;
    }
    if (h_joined != nboxes) {
        amrex::Print() << "FAIL nested " << name << ": work queued after the loop saw "
                       << h_joined << " of " << nboxes << " loop results\n";
        ++nfail;
    }
    if (nfail == 0) {
        amrex::Print() << "PASS nested " << name << " (origin stream " << origin << ")\n";
    }
    return nfail;
}

enum struct NoForkKind { MFIterHost, MFIterAllBoxes, StreamIterDeviceSync };
constexpr NoForkKind all_nofork_kinds[] = {
    NoForkKind::MFIterHost, NoForkKind::MFIterAllBoxes, NoForkKind::StreamIterDeviceSync
};

const char* nofork_name (NoForkKind kind)
{
    switch (kind) {
    case NoForkKind::MFIterHost: return "MFIter outside launch region";
    case NoForkKind::MFIterAllBoxes: return "MFIter(AllBoxes)";
    case NoForkKind::StreamIterDeviceSync: return "StreamIter(device_sync)";
    }
    return "";
}

// Runs a loop that does not fork onto other streams from inside an MFIter
// iteration on stream 1, between a slow producer and a consumer queued on
// that stream, and checks that the consumer sees the producer's result.
int test_nested_nofork_loop (MultiFab const& mf, Long iters_per_ms, NoForkKind kind,
                             bool warmup = false)
{
    const char* name = nofork_name(kind);
    const int nboxes = mf.local_size();

    Gpu::DeviceVector<int> written(1, 0);
    Gpu::DeviceVector<int> seen(1, -1);
    Gpu::DeviceVector<unsigned> sink(1, 0);
    auto* pw = written.data();
    auto* ps = seen.data();
    auto* pk = sink.data();
    const Long long_wait = warmup ? 0 : 50 * iters_per_ms;
    int origin = -1;
    int after = -1;
    int count = 0;
    Gpu::synchronize();

    {
        Gpu::StreamOrderedRegion stream_ordered;
        for (MFIter omfi(mf); omfi.isValid(); ++omfi) {
            if (omfi.LocalIndex() != 1) { continue; }
            origin = Gpu::Device::streamIndex();

            amrex::single_task([=] AMREX_GPU_DEVICE () {
                pk[0] = spin(long_wait);
                pw[0] = 1;
            });

            if (kind == NoForkKind::MFIterHost) {
                Gpu::LaunchSafeGuard lsg(false);
                for (MFIter mfi(mf); mfi.isValid(); ++mfi) { ++count; }
            } else if (kind == NoForkKind::MFIterAllBoxes) {
                const unsigned char flags = MFIter::AllBoxes;
                for (MFIter mfi(mf, flags); mfi.isValid(); ++mfi) { ++count; }
            } else {
                Gpu::StreamItInfo info;
                info.device_sync = true;
                for (Gpu::StreamIter sit(nboxes, info); sit.isValid(); ++sit) { ++count; }
            }

            after = Gpu::Device::streamIndex();
            amrex::single_task([=] AMREX_GPU_DEVICE () { ps[0] = pw[0]; });
        }
    }
    Gpu::synchronize();
    if (warmup) { return 0; }

    int h_seen = 0;
    Gpu::copy(Gpu::deviceToHost, seen.begin(), seen.end(), &h_seen);

    int nfail = 0;
    if (after != origin) {
        amrex::Print() << "FAIL nested " << name << ": current stream after the loop is "
                       << after << ", expected " << origin << "\n";
        ++nfail;
    }
    if (h_seen != 1) {
        amrex::Print() << "FAIL nested " << name << ": work queued after the loop ran "
                       << "before the work queued ahead of it on stream " << origin << "\n";
        ++nfail;
    }
    if (nfail == 0) {
        amrex::Print() << "PASS nested " << name << " (origin stream " << origin
                       << ", " << count << " items)\n";
    }
    return nfail;
}

// Builds the arrays of a new MultiFab in iteration 1 of a StreamIter loop (or
// an MFIter loop, after a nested MFIter has finished in that iteration),
// behind a slow kernel on stream 1, and reads them in iteration 2.
int test_loop_arrays (BoxArray const& ba, Long iters_per_ms, bool outer_mfiter,
                      bool warmup = false)
{
    MultiFab mf(ba, DistributionMapping(ba), 1, 0);
    const int n = mf.local_size();
    Gpu::DeviceVector<int> lo(n, -1);
    Gpu::DeviceVector<unsigned> sink(1, 0);
    auto* plo = lo.data();
    auto* pk = sink.data();
    const Long long_wait = warmup ? 0 : 50 * iters_per_ms;
    Gpu::synchronize();

    auto body = [&] (int i) {
        if (i == 1) {
            if (outer_mfiter) {
                // A nested loop that has finished must not hide the outer one.
                for (MFIter mfi(mf); mfi.isValid(); ++mfi) {}
            }
            amrex::single_task([=] AMREX_GPU_DEVICE () { pk[0] = spin(long_wait); });
            amrex::ignore_unused(mf.const_arrays());  // first use: copied to the device
        } else if (i == 2) {
            auto const& ma = mf.const_arrays();
            amrex::single_task([=] AMREX_GPU_DEVICE () {
                for (int b = 0; b < n; ++b) { plo[b] = ma[b].begin[0]; }
            });
        }
    };
    {
        Gpu::StreamOrderedRegion stream_ordered;
        if (outer_mfiter) {
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) { body(mfi.LocalIndex()); }
        } else {
            for (Gpu::StreamIter sit(n); sit.isValid(); ++sit) { body(sit()); }
        }
    }
    Gpu::synchronize();
    if (warmup) { return 0; }
    const char* name = outer_mfiter ? "MFIter arrays after nested MFIter" : "StreamIter arrays";

    Vector<int> h_lo(n);
    Gpu::copy(Gpu::deviceToHost, lo.begin(), lo.end(), h_lo.begin());
    int nwrong = 0;
    for (int b = 0; b < n; ++b) {
        if (h_lo[b] != mf.fabbox(mf.IndexArray()[b]).smallEnd(0)) { ++nwrong; }
    }
    if (nwrong > 0) {
        amrex::Print() << "FAIL " << name << ": " << nwrong << " of " << n
                       << " Array4s read before their copy to the device\n";
        return 1;
    }
    amrex::Print() << "PASS " << name << "\n";
    return 0;
}

// Writes a MultiFab with VisMF::AsyncWrite (amrex.async_out=1) in a
// StreamOrderedRegion right after filling it, then reads the file back. A
// smoke test: AsyncWrite's own min/max reduction waits for the device before
// its device-to-host copies, and those copies usually finish before the
// writer thread reads them, so a missing wait in AsyncOut::Submit is not
// reliably caught here.
int test_async_write (BoxArray const& ba)
{
    if (!AsyncOut::UseAsyncOut()) {
        amrex::Print() << "SKIP async write: amrex.async_out is off\n";
        return 0;
    }
    const std::string dir = "StreamOrderedRegion_async";
    const std::string name = dir + "/mf";
    amrex::UtilCreateCleanDirectory(dir, true);
    constexpr Real value = 42.0;
    MultiFab mf(ba, DistributionMapping(ba), 1, 0);
    mf.setVal(0.0);
    Gpu::synchronize();

    {
        Gpu::StreamOrderedRegion stream_ordered;
        mf.setVal(value);
        VisMF::AsyncWrite(mf, name);
    }
    AsyncOut::Finish();

    MultiFab in;
    VisMF::Read(in, name);
    const Real lo = in.min(0);
    const Real hi = in.max(0);
    ParallelDescriptor::Barrier();
    if (ParallelDescriptor::IOProcessor()) { FileSystem::RemoveAll(dir); }
    if (lo != value || hi != value) {
        amrex::Print() << "FAIL async write: file holds values in [" << lo << ", " << hi
                       << "], expected " << value << "\n";
        return 1;
    }
    amrex::Print() << "PASS async write\n";
    return 0;
}

#ifdef AMREX_USE_OMP
// An outer MFIter run by an active OpenMP team outside the launch region, with
// an MFIter nested in each iteration. With 5 boxes and 2 threads, the threads
// run 3 and 2 outer iterations. A worksharing construct in the nested MFIter
// would then be reached a different number of times by each thread and hang.
int test_nested_host_mfiter_omp ()
{
    BoxArray ba(Box(IntVect(0), IntVect(AMREX_D_DECL(19,3,3))));
    ba.maxSize(4);  // 5 boxes
    MultiFab outer(ba, DistributionMapping(ba), 1, 0);
    MultiFab inner(ba, DistributionMapping(ba), 1, 0);
    long count = 0;
    for (int rep = 0; rep < 10; ++rep) {
        Gpu::StreamOrderedRegion stream_ordered;
        Gpu::LaunchSafeGuard lsg(false);
#pragma omp parallel num_threads(2) reduction(+:count)
        for (MFIter omfi(outer); omfi.isValid(); ++omfi) {
            for (MFIter mfi(inner); mfi.isValid(); ++mfi) { ++count; }
        }
    }
    amrex::Print() << "PASS nested host MFIter in OpenMP team (" << count << " items)\n";
    return 0;
}
#endif

#ifdef AMREX_USE_CUDA
// Captures, into a CUDA graph, an MFIter over 2 boxes whose iterations build
// and read FabArray::arrays() for the first time. Only the loop's streams may
// join the capture: a stream the loop does not use would stay in it until the
// capture ends, so that work queued on it meanwhile, e.g., by another thread,
// would be captured too or fail. (CUDA does not report such a stream as
// unjoined when no work was queued on it.)
//
// With synced_mf, a loop over its boxes that synchronizes the host
// (MFItInfo::SetDeviceSync(true)) runs first in the same region: it must not
// leave its streams counted as in use after it ends. With all_boxes, the
// captured loop is an MFIter(AllBoxes), which uses its streams like any loop.
int test_graph_capture_arrays (MultiFab const* synced_mf = nullptr, bool all_boxes = false,
                               bool warmup = false)
{
    BoxArray ba(Box(IntVect(1000), IntVect(AMREX_D_DECL(1007,1003,1003))));
    ba.maxSize(4);  // 2 boxes
    MultiFab mf(ba, DistributionMapping(ba), 1, 0);
    const int n = mf.local_size();
    if (n >= Gpu::Device::streamPoolSize()) {
        amrex::Print() << "SKIP graph capture: needs more GPU streams than local boxes\n";
        return 0;
    }
    Gpu::DeviceVector<int> lo(n, -1);
    auto* plo = lo.data();
    Gpu::synchronize();

    cudaGraph_t graph = nullptr;
    cudaError_t end_status = cudaSuccess;
    int ncaptured = 0;  // streams the loop does not use found in the capture
    {
        Gpu::StreamOrderedRegion stream_ordered;
        if (synced_mf) {
            for (MFIter mfi(*synced_mf, MFItInfo().SetDeviceSync(true)); mfi.isValid(); ++mfi) {}
        }
        const cudaStream_t s = Gpu::gpuStream();
        if (!warmup) {
            AMREX_CUDA_SAFE_CALL(cudaStreamBeginCapture(s, cudaStreamCaptureModeRelaxed));
        }
        const int nused = n;  // streams the loop uses
        auto body = [&] (int b) {
            auto const& ma = mf.const_arrays();
            amrex::single_task([=] AMREX_GPU_DEVICE () { plo[b] = ma[b].begin[0]; });
            if (!warmup) {
                const int idx = Gpu::Device::streamIndex();
                for (int i = nused; i < Gpu::Device::streamPoolSize(); ++i) {
                    Gpu::Device::setStreamIndex(i);
                    cudaStreamCaptureStatus status{};
                    AMREX_CUDA_SAFE_CALL(cudaStreamIsCapturing(Gpu::gpuStream(), &status));
                    if (status != cudaStreamCaptureStatusNone) { ++ncaptured; }
                }
                Gpu::Device::setStreamIndex(idx);
            }
        };
        if (all_boxes) {
            // One rank: the global index is the local one.
            const unsigned char flags = MFIter::AllBoxes;
            for (MFIter mfi(mf, flags); mfi.isValid(); ++mfi) { body(mfi.index()); }
        } else {
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) { body(mfi.LocalIndex()); }
        }
        if (!warmup) {
            end_status = cudaStreamEndCapture(s, &graph);
            (void)cudaGetLastError();
        }
    }
    if (warmup) {
        Gpu::synchronize();
        return 0;
    }
    if (end_status != cudaSuccess) {
        amrex::Print() << "FAIL graph capture: ending the capture failed with "
                       << cudaGetErrorName(end_status) << "\n";
        if (graph) { cudaGraphDestroy(graph); }
        return 1;
    }
    if (ncaptured > 0) {
        amrex::Print() << "FAIL graph capture: streams the loop does not use joined the "
                       << "capture (" << ncaptured << " times)\n";
        if (graph) { cudaGraphDestroy(graph); }
        return 1;
    }

    cudaGraphExec_t exec = nullptr;
    AMREX_CUDA_SAFE_CALL(cudaGraphInstantiate(&exec, graph, 0));
    AMREX_CUDA_SAFE_CALL(cudaGraphLaunch(exec, Gpu::gpuStream()));
    Gpu::synchronize();
    AMREX_CUDA_SAFE_CALL(cudaGraphExecDestroy(exec));
    AMREX_CUDA_SAFE_CALL(cudaGraphDestroy(graph));

    Vector<int> h_lo(n);
    Gpu::copy(Gpu::deviceToHost, lo.begin(), lo.end(), h_lo.begin());
    int nwrong = 0;
    for (int b = 0; b < n; ++b) {
        if (h_lo[b] != mf.fabbox(mf.IndexArray()[b]).smallEnd(0)) { ++nwrong; }
    }
    if (nwrong > 0) {
        amrex::Print() << "FAIL graph capture: " << nwrong << " of " << n
                       << " Array4s read wrong in the graph\n";
        return 1;
    }
    amrex::Print() << "PASS graph capture of arrays in a 2-box "
                   << (all_boxes ? "MFIter(AllBoxes) loop" : "loop")
                   << (synced_mf ? " after a synchronizing loop" : "") << "\n";
    return 0;
}
#endif

// An outer MFIter outside any Gpu::StreamOrderedRegion (it synchronizes, and
// runs its iterations on several streams). In iteration 1, a region with a
// one-stream MFIter builds the arrays of a new MultiFab behind a slow kernel;
// iteration 2, on another stream, reads them. The arrays must be ordered for
// the outer loop's streams, not just the inner loop's one.
int test_arrays_in_region_in_synced_loop (BoxArray const& ba, Long iters_per_ms,
                                          bool warmup = false)
{
    MultiFab outer(ba, DistributionMapping(ba), 1, 0);
    BoxArray ba2 = ba;
    if (!warmup) { ba2.shift(IntVect(2000)); }  // differ from the warm-up's
    MultiFab mf(ba2, DistributionMapping(ba2), 1, 0);
    const int n = mf.local_size();
    Gpu::DeviceVector<int> lo(n, -1);
    Gpu::DeviceVector<unsigned> sink(1, 0);
    auto* plo = lo.data();
    auto* pk = sink.data();
    const Long long_wait = warmup ? 0 : 50 * iters_per_ms;
    int stream1 = -1;
    int stream2 = -1;
    Gpu::synchronize();

    for (MFIter omfi(outer); omfi.isValid(); ++omfi) {
        const int i = omfi.LocalIndex();
        if (i == 1) {
            stream1 = Gpu::Device::streamIndex();
            Gpu::StreamOrderedRegion stream_ordered;
            for (MFIter mfi(mf, MFItInfo().SetNumStreams(1)); mfi.isValid(); ++mfi) {
                if (mfi.LocalIndex() != 0) { continue; }
                amrex::single_task([=] AMREX_GPU_DEVICE () { pk[0] = spin(long_wait); });
                amrex::ignore_unused(mf.const_arrays());  // first use: copied to the device
            }
        } else if (i == 2) {
            stream2 = Gpu::Device::streamIndex();
            auto const& ma = mf.const_arrays();
            amrex::single_task([=] AMREX_GPU_DEVICE () {
                for (int b = 0; b < n; ++b) { plo[b] = ma[b].begin[0]; }
            });
        }
    }
    Gpu::synchronize();
    if (warmup) { return 0; }

    Vector<int> h_lo(n);
    Gpu::copy(Gpu::deviceToHost, lo.begin(), lo.end(), h_lo.begin());
    int nwrong = 0;
    for (int b = 0; b < n; ++b) {
        if (h_lo[b] != mf.fabbox(mf.IndexArray()[b]).smallEnd(0)) { ++nwrong; }
    }
    if (nwrong > 0) {
        amrex::Print() << "FAIL arrays in region in synchronizing loop: " << nwrong << " of "
                       << n << " Array4s read on stream " << stream2
                       << " before their copy to the device\n";
        return 1;
    }
    amrex::Print() << "PASS arrays in region in synchronizing loop (streams " << stream1
                   << ", " << stream2 << ")\n";
    return 0;
}

// In a Gpu::StreamOrderedRegion with Gpu::setNoSyncRegion(false), a plain
// StreamIter must synchronize at its end like outside the region: the host
// reads the results of device-to-host copies queued in the loop right after
// it, without waiting.
int test_streamiter_sync_turned_on (Long iters_per_ms, bool warmup = false)
{
    constexpr int n = 4;
    Gpu::DeviceVector<int> d(n, 0);
    Gpu::DeviceVector<unsigned> sink(n, 0);
    auto* pd = d.data();
    auto* pk = sink.data();
    auto* h = static_cast<int*>(The_Pinned_Arena()->alloc(n * sizeof(int)));
    for (int i = 0; i < n; ++i) { h[i] = 0; }
    const Long long_wait = warmup ? 0 : 50 * iters_per_ms;
    Gpu::synchronize();

    Vector<int> seen(n, -1);
    {
        Gpu::StreamOrderedRegion stream_ordered;
        const bool prev = Gpu::setNoSyncRegion(false);
        for (Gpu::StreamIter sit(n); sit.isValid(); ++sit) {
            const int i = sit();
            amrex::single_task([=] AMREX_GPU_DEVICE () {
                pk[i] = spin(long_wait);
                pd[i] = i + 1;
            });
            Gpu::dtoh_memcpy_async(h + i, pd + i, sizeof(int));
        }
        for (int i = 0; i < n; ++i) { seen[i] = h[i]; }  // no explicit wait
        Gpu::setNoSyncRegion(prev);
    }
    Gpu::synchronize();
    The_Pinned_Arena()->free(h);
    if (warmup) { return 0; }

    int nstale = 0;
    for (int i = 0; i < n; ++i) { if (seen[i] != i + 1) { ++nstale; } }
    if (nstale > 0) {
        amrex::Print() << "FAIL StreamIter with no-sync turned off: " << nstale << " of " << n
                       << " results read before the loop's copies finished\n";
        return 1;
    }
    amrex::Print() << "PASS StreamIter with no-sync turned off\n";
    return 0;
}

// An MFIter with no local items, constructed in an iteration on stream 1 of a
// Gpu::StreamOrderedRegion, must keep stream 1 current: work queued before
// its iterations (e.g., the fill of a FillPatchIterator) would otherwise go
// to a stream that nothing orders with the work around the loop.
int test_empty_loop_keeps_stream (MultiFab const& mf)
{
    MultiFab empty(BoxArray{}, DistributionMapping{}, 1, 0);
    int origin = -1;
    int during = -1;
    {
        Gpu::StreamOrderedRegion stream_ordered;
        for (MFIter omfi(mf); omfi.isValid(); ++omfi) {
            if (omfi.LocalIndex() != 1) { continue; }
            origin = Gpu::Device::streamIndex();
            MFIter mfi(empty);
            during = Gpu::Device::streamIndex();
        }
    }
    Gpu::synchronize();
    if (during != origin) {
        amrex::Print() << "FAIL empty MFIter: current stream changed from " << origin
                       << " to " << during << "\n";
        return 1;
    }
    amrex::Print() << "PASS empty MFIter keeps the current stream (" << origin << ")\n";
    return 0;
}

// Forwards to The_Arena and counts frees of watched memory whose kernel has
// not finished yet.
class WatchingArena final : public Arena
{
public:
    void* alloc (std::size_t sz) override { return The_Arena()->alloc(sz); }

    void free (void* p) override {
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            auto it = m_done.find(p);
            if (it != m_done.end()) {
                if (*(it->second) == 0) { ++m_early; }
                m_done.erase(it);
                ++m_freed;
            }
        }
        The_Arena()->free(p);
    }

    void watch (void* p, volatile int* done) {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_done[p] = done;
    }

    int early () {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_early;
    }

    int freed () {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_freed;
    }

private:
    std::mutex m_mutex;
    std::map<void*, volatile int*> m_done;
    int m_early = 0;
    int m_freed = 0;
};

// One thread queues kernels on stream 1 and releases their memory with
// Gpu::freeAsync while the main thread joins stream 1 into stream 0.
// Synchronizing stream 0 afterwards must only free memory whose kernel is done.
#ifdef AMREX_USE_OMP
int test_join_concurrent_free (Long iters_per_ms, int ntrials)
{
    constexpr int nbufs = 32;
    WatchingArena arena;
    auto* done = static_cast<int*>(The_Pinned_Arena()->alloc(nbufs * sizeof(int)));
    const Long wait = iters_per_ms;
    int nmid = 0;  // trials in which the join split the worker's frees

    for (int trial = 0; trial < ntrials; ++trial) {
        std::vector<void*> bufs(nbufs);
        for (int k = 0; k < nbufs; ++k) {
            done[k] = 0;
            bufs[k] = arena.alloc(sizeof(unsigned));
            arena.watch(bufs[k], done + k);
        }
        Gpu::synchronize();

        // OpenMP threads, not a std::thread: AMReX keeps the current stream
        // per OpenMP thread, so a std::thread would share the main thread's.
        std::atomic<bool> ready{false};
        std::atomic<bool> go{false};
#pragma omp parallel num_threads(2)
        if (omp_get_thread_num() == 1) {
            ready = true;
            while (!go) {}
            Gpu::Device::setStreamIndex(1);
            for (int k = 0; k < nbufs; ++k) {
                auto* p = static_cast<unsigned*>(bufs[k]);
                volatile int* d = done + k;
                amrex::single_task([=] AMREX_GPU_DEVICE () {
                    *p = spin(wait);
                    *d = 1;
                });
                Gpu::freeAsync(&arena, bufs[k]);
            }
            Gpu::Device::setStreamIndex(0);
        } else {
            while (!ready) {}
            go = true;
            // Join at a different point of the worker's queueing in each trial.
            const auto t0 = std::chrono::steady_clock::now();
            const auto delay = std::chrono::microseconds(5 * (trial % 40));
            while (std::chrono::steady_clock::now() - t0 < delay) {}
            Gpu::Fan fan;
            Gpu::theStreamSequencer().fork(fan, {.nitems = 2, .max_width = 2,
                                           .fork = Gpu::StreamEdge::Gpu,
                                           .join = Gpu::StreamEdge::Gpu});
            Gpu::theStreamSequencer().join(fan);
        }

        Gpu::Device::setStreamIndex(0);
        const int freed0 = arena.freed();
        Gpu::streamSynchronize();  // frees what the join handed to stream 0
        const int moved = arena.freed() - freed0;
        if (moved > 0 && moved < nbufs) { ++nmid; }
        Gpu::synchronize();
    }

    The_Pinned_Arena()->free(done);

    const int early = arena.early();
    if (early > 0) {
        amrex::Print() << "FAIL join with concurrent freeAsync: " << early
                       << " buffers freed while their kernel was still running ("
                       << ntrials << " trials, " << nmid << " joined mid-burst)\n";
        return 1;
    }
    amrex::Print() << "PASS join with concurrent freeAsync (" << ntrials << " trials, "
                   << nmid << " joined mid-burst)\n";
    return 0;
}
#endif

}

int main (int argc, char* argv[])
{
    int nfail = 0;
    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, [] () {
        ParmParse pp("amrex");
        int async_out = 1;  // for test_async_write
        pp.queryAdd("async_out", async_out);
    });
    {
        if (Gpu::Device::numGpuStreams() < 2) {
            amrex::Print() << "SKIP: needs at least 2 GPU streams\n";
        } else {
            int ntrials = 100;
            ParmParse pp;
            pp.query("ntrials", ntrials);

            const Long iters_per_ms = spin_iters_per_ms();

            MFIter::allowMultipleMFIters(true);
            BoxArray ba(Box(IntVect(0), IntVect(15)));
            ba.maxSize(4);  // at least four boxes in any dimension
            MultiFab mf(ba, DistributionMapping(ba), 1, 0);

            if (mf.local_size() < 2) {
                amrex::Print() << "SKIP nested loops: needs at least 2 local boxes\n";
            } else {
                // Launch every kernel once first: loading a kernel's module on its
                // first launch (e.g., CUDA lazy loading) can wait for the device.
                for (auto kind : all_loop_kinds) {
                    test_nested_loop(mf, iters_per_ms, kind, true);
                }
                for (auto kind : all_nofork_kinds) {
                    test_nested_nofork_loop(mf, iters_per_ms, kind, true);
                }
                for (auto kind : all_loop_kinds) {
                    nfail += test_nested_loop(mf, iters_per_ms, kind);
                }
                for (auto kind : all_nofork_kinds) {
                    nfail += test_nested_nofork_loop(mf, iters_per_ms, kind);
                }
                nfail += test_empty_loop_keeps_stream(mf);
            }

            if (mf.local_size() < 3) {
                amrex::Print() << "SKIP loop arrays: needs at least 3 local boxes\n";
            } else {
                // The boxes differ from the warm-up's, so that a device table
                // left over from it cannot pass for the new one.
                for (bool outer_mfiter : {false, true}) {
                    test_loop_arrays(ba, iters_per_ms, outer_mfiter, true);
                    BoxArray ba2 = ba;
                    ba2.shift(IntVect(1000));
                    nfail += test_loop_arrays(ba2, iters_per_ms, outer_mfiter);
                }
            }
            nfail += test_async_write(ba);
#ifdef AMREX_USE_OMP
            if (OpenMP::get_max_threads() < 2) {
                amrex::Print() << "SKIP join with concurrent freeAsync: needs 2 OpenMP threads\n";
            } else {
                nfail += test_join_concurrent_free(iters_per_ms, ntrials);
            }
#else
            amrex::Print() << "SKIP join with concurrent freeAsync: needs OpenMP\n";
#endif
#ifdef AMREX_USE_OMP
            nfail += test_nested_host_mfiter_omp();
#endif
            // Before the graph capture tests: a loop-stream count leaked by
            // one of them would hide a failure here.
            if (mf.local_size() < 3) {
                amrex::Print() << "SKIP arrays in region in synchronizing loop: "
                               << "needs at least 3 local boxes\n";
            } else {
                test_arrays_in_region_in_synced_loop(ba, iters_per_ms, true);
                nfail += test_arrays_in_region_in_synced_loop(ba, iters_per_ms);
            }
            test_streamiter_sync_turned_on(iters_per_ms, true);
            nfail += test_streamiter_sync_turned_on(iters_per_ms);
#ifdef AMREX_USE_CUDA
            test_graph_capture_arrays(nullptr, false, true);
            nfail += test_graph_capture_arrays();
            nfail += test_graph_capture_arrays(&mf);
            if (ParallelDescriptor::NProcs() == 1) {
                nfail += test_graph_capture_arrays(nullptr, true);
            }
#endif
        }
    }
    amrex::Finalize();
    return nfail == 0 ? 0 : 1;
}
