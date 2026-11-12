#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>
#include <AMReX_Loom.H>
#include <AMReX_Print.H>

#include <cstdio>
#include <vector>

using namespace amrex;

/**
 * LOOM backend: each translation unit that launches kernels embeds its own
 * code object (libamrex's units and this file), and the runtime loads one
 * only at the first launch of one of its kernels, never at startup.
 *
 * The kernels also take by-value arguments narrower than a kernel argument
 * slot (char, short, bool, unsigned char) mixed with wider ones, which must
 * reach the device exactly.
 */

namespace {
int failures = 0;
void report (bool ok, const char* what)
{
    std::printf("%s: %s\n", ok ? "ok" : "FAIL", what);
    if (!ok) { ++failures; }
}
}

int main (int argc, char* argv[])
{
    // The runtime is up and nothing has launched: no code object is loaded,
    // although the program embeds dozens.
    report(amrex_loom::device_count() == 1, "one device");
    report(amrex_loom::loaded_code_object_count() == 0,
           "no code object is loaded before the first launch");

    amrex::Initialize(argc, argv);
    {
        const std::size_t at_init = amrex_loom::loaded_code_object_count();
        amrex::Print() << "code objects loaded by amrex::Initialize: " << at_init << "\n";

        constexpr int n = 256;
        auto* d = static_cast<long long*>(The_Arena()->alloc(n * sizeof(long long)));
        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) { d[i] = i; });
        Gpu::streamSynchronize();
        const std::size_t after_first = amrex_loom::loaded_code_object_count();
        report(after_first == at_init + 1,
               "the first launch from this file loads this file's code object");

        // Not const, so the lambda captures them instead of their values.
        char c = -3;
        short s = 1234;
        bool b = true;
        unsigned char u = 200;
        double x = 0.5;
        char c2 = 7;
        int k = -100000;
        ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
        {
            d[i] = d[i] * 1000000 + c + s + (b ? 10 : 0) + u + static_cast<long long>(4 * x)
                 + c2 + k;
        });
        Gpu::streamSynchronize();
        report(amrex_loom::loaded_code_object_count() == after_first,
               "another kernel of this file reuses the loaded code object");

        std::vector<long long> h(n);
        Gpu::dtoh_memcpy(h.data(), d, n * sizeof(long long));
        bool values = true;
        for (int i = 0; i < n; ++i) {
            const long long expected = i * 1000000LL + c + s + 10 + u + 2 + c2 + k;
            if (h[i] != expected) {
                if (values) {
                    std::printf("d[%d] = %lld, expected %lld\n", i, h[i], expected);
                }
                values = false;
            }
        }
        report(values, "sub-word and wider by-value arguments reach the kernel");
        The_Arena()->free(d);
    }
    amrex::Finalize();

    std::printf("%s\n", failures ? "FAILED" : "PASS");
    return failures ? 1 : 0;
}
