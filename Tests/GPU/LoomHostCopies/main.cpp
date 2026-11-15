#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>

#include <cstdio>
#include <vector>

using namespace amrex;

/**
 * LOOM backend: asynchronous copies between device memory and pageable host
 * memory behave as hipMemcpyAsync makes them. HRX reads and writes the host
 * range when a transfer executes, and the runtime batches stream work, so
 * without the runtime's staging a reused source or an early read would see
 * the wrong data:
 *
 * - host to device: the source may be overwritten as soon as the call returns;
 * - device to host: the destination holds the data when the call returns;
 * - pinned host memory stays asynchronous and is complete after a stream
 *   synchronization.
 */

namespace {
int failures = 0;
void report (bool ok, const char* what)
{
    std::printf("%s: %s\n", ok ? "ok" : "FAIL", what);
    if (!ok) { ++failures; }
}

template <class Check>
bool all_of (int n, Check check)
{
    for (int i = 0; i < n; ++i) {
        if (!check(i)) {
            std::printf("first mismatch at %d\n", i);
            return false;
        }
    }
    return true;
}
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        constexpr int n = 1 << 16;
        constexpr std::size_t bytes = n * sizeof(long long);
        auto* d = static_cast<long long*>(The_Arena()->alloc(bytes));

        // Host to device from pageable memory, source reused at once.
        {
            std::vector<long long> h(n);
            for (int i = 0; i < n; ++i) { h[i] = 3LL * i + 1; }
            Gpu::htod_memcpy_async(d, h.data(), bytes);
            for (auto& v : h) { v = -1; }
            h.clear();
            h.shrink_to_fit();
            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) { d[i] += 1; });
            std::vector<long long> back(n);
            Gpu::dtoh_memcpy(back.data(), d, bytes);
            report(all_of(n, [&] (int i) { return back[i] == 3LL * i + 2; }),
                   "a pageable host-to-device source may be reused when the call returns");
        }

        // Device to host into pageable memory, read without a synchronization.
        {
            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) { d[i] = 7LL * i; });
            std::vector<long long> h(n, -1);
            Gpu::dtoh_memcpy_async(h.data(), d, bytes);
            report(all_of(n, [&] (int i) { return h[i] == 7LL * i; }),
                   "a pageable device-to-host destination is complete when the call returns");
        }

        // Pinned memory: asynchronous, complete after a stream synchronization.
        {
            auto* p = static_cast<long long*>(The_Pinned_Arena()->alloc(bytes));
            for (int i = 0; i < n; ++i) { p[i] = 5LL * i; }
            Gpu::htod_memcpy_async(d, p, bytes);
            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) { d[i] *= 2; });
            Gpu::dtoh_memcpy_async(p, d, bytes);
            Gpu::streamSynchronize();
            report(all_of(n, [&] (int i) { return p[i] == 10LL * i; }),
                   "pinned copies are complete after a stream synchronization");
            The_Pinned_Arena()->free(p);
        }

        The_Arena()->free(d);
    }
    amrex::Finalize();

    std::printf("%s\n", failures == 0 ? "PASSED" : "FAILED");
    return failures == 0 ? 0 : 1;
}
