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
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <atomic>
#include <chrono>
#include <map>
#include <mutex>
#include <thread>
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
 * 4. Memory released with Gpu::freeAsync on a stream while another thread
 *    joins that stream (Gpu::Device::joinStreams) must not be freed before
 *    the kernels queued ahead of it on that stream finish. This is a stress
 *    test: the window it probes is narrow, and on CUDA it is usually
 *    shorter than a kernel launch, so it mainly catches regressions on HIP.
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

        std::atomic<bool> ready{false};
        std::atomic<bool> go{false};
        std::thread worker([&] () {
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
        });
        while (!ready) {}
        go = true;
        // Join at a different point of the worker's queueing in each trial.
        const auto t0 = std::chrono::steady_clock::now();
        const auto delay = std::chrono::microseconds(5 * (trial % 40));
        while (std::chrono::steady_clock::now() - t0 < delay) {}
        Gpu::Device::joinStreams(2, 0);
        worker.join();

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
            nfail += test_join_concurrent_free(iters_per_ms, ntrials);
        }
    }
    amrex::Finalize();
    return nfail == 0 ? 0 : 1;
}
