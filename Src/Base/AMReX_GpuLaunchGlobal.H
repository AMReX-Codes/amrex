#ifndef AMREX_GPU_LAUNCH_GLOBAL_H_
#define AMREX_GPU_LAUNCH_GLOBAL_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>

#if defined(AMREX_USE_CUDA) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

namespace amrex {

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    // We cannot take rvalue lambdas.
    template<class L>
    AMREX_GPU_GLOBAL void launch_global (L f0) { f0(); }

#if defined(AMREX_USE_CUDA) && defined(AMREX_GPU_MIN_BLOCKS)
    namespace detail {
        // Min blocks of mt threads for AMREX_GPU_MIN_BLOCKS blocks of
        // AMREX_GPU_MAX_THREADS threads. Rounded down, so for mt up to
        // AMREX_GPU_MIN_BLOCKS*AMREX_GPU_MAX_THREADS the kernel's register cap is
        // not below -maxrregcount, which non-inlined callees are compiled with.
        // At most 16, the lowest per-SM block limit of supported GPUs.
        constexpr int launch_min_blocks (int mt) {
            int n = AMREX_GPU_MIN_BLOCKS * AMREX_GPU_MAX_THREADS / mt;
            return (n < 1) ? 1 : ((n > 16) ? 16 : n);
        }
    }
#endif

    // With AMREX_GPU_MIN_BLOCKS (CUDA only), registers per thread are capped so
    // that this many blocks of AMREX_GPU_MAX_THREADS threads fit on an SM.
    // Note that -maxrregcount is ignored for kernels with launch bounds.
    template<int amrex_launch_bounds_max_threads, class L>
#if defined(AMREX_USE_CUDA) && defined(AMREX_GPU_MIN_BLOCKS)
    __launch_bounds__(amrex_launch_bounds_max_threads,
                      detail::launch_min_blocks(amrex_launch_bounds_max_threads))
#else
    __launch_bounds__(amrex_launch_bounds_max_threads)
#endif
    AMREX_GPU_GLOBAL void launch_global (L f0) { f0(); }

    template<int amrex_launch_bounds_max_threads, int min_blocks, class L>
    __launch_bounds__(amrex_launch_bounds_max_threads, min_blocks)
    AMREX_GPU_GLOBAL void launch_global (L f0) { f0(); }
#endif

}

#if defined(AMREX_USE_CUDA) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#endif
