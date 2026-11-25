#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>

#include <cstdio>

using namespace amrex;

/**
 * LOOM backend: pinned host memory (The_Pinned_Arena) that kernels read and
 * write directly stays coherent across AMReX's completion points.
 *
 * HRX maps pinned buffers into the GPU's L2 cache. A kernel that reads one
 * leaves its lines in L2; if the host then rewrites the buffer, a later
 * kernel must not see the old lines. AMReX only lets the host touch pinned
 * memory after a synchronization, so the runtime must make every
 * synchronization invalidate (and write back) L2 after the last kernel,
 * whatever submission mode it uses.
 *
 * Each round, the host rewrites a small pinned input, a kernel copies it into
 * device memory, a second kernel writes a pinned output, and the host checks
 * the output after a synchronization. A stale L2 line shows up as last
 * round's value.
 */

namespace {
int failures = 0;
void report (bool ok, const char* what)
{
    std::printf("%s: %s\n", ok ? "ok" : "FAIL", what);
    if (!ok) { ++failures; }
}

// Returns the number of rounds whose output did not match.
int rounds (int nrounds, int n, int* in, int* dev, int* out, bool all_streams)
{
    int bad = 0;
    for (int r = 0; r < nrounds; ++r) {
        for (int i = 0; i < n; ++i) { in[i] = r * n + i; }
        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) { dev[i] = in[i]; });
        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) { out[i] = dev[i] + 1; });
        if (all_streams) {
            Gpu::Device::synchronize();
        } else {
            Gpu::streamSynchronize();
        }
        for (int i = 0; i < n; ++i) {
            if (out[i] != r * n + i + 1) {
                if (bad == 0) {
                    std::printf("round %d: out[%d] = %d, expected %d\n", r, i, out[i],
                                r * n + i + 1);
                }
                ++bad;
                break;
            }
        }
    }
    return bad;
}
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        // Small enough to stay resident in L2 between rounds.
        constexpr int n = 4096;
        constexpr int nrounds = 500;
        auto* in = static_cast<int*>(The_Pinned_Arena()->alloc(n * sizeof(int)));
        auto* out = static_cast<int*>(The_Pinned_Arena()->alloc(n * sizeof(int)));
        auto* dev = static_cast<int*>(The_Arena()->alloc(n * sizeof(int)));

        int bad = rounds(nrounds, n, in, dev, out, false);
        if (bad) { std::printf("%d of %d rounds read stale data\n", bad, nrounds); }
        report(bad == 0,
               "kernels see host writes to pinned memory after a stream synchronization");

        // The same on another stream, synchronized through the whole device.
        if (Gpu::numGpuStreams() > 1) {
            Gpu::Device::setStreamIndex(1);
            bad = rounds(nrounds, n, in, dev, out, true);
            Gpu::Device::resetStreamIndex();
            if (bad) { std::printf("%d of %d rounds read stale data\n", bad, nrounds); }
            report(bad == 0,
                   "the same on a second stream after a device synchronization");
        }

        The_Arena()->free(dev);
        The_Pinned_Arena()->free(out);
        The_Pinned_Arena()->free(in);
    }
    amrex::Finalize();

    std::printf("%s\n", failures == 0 ? "PASSED" : "FAILED");
    return failures == 0 ? 0 : 1;
}
