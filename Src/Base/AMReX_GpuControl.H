#ifndef AMREX_GPU_CONTROL_H_
#define AMREX_GPU_CONTROL_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuTypes.H>

#include <utility>

#if defined(AMREX_USE_CUDA) && (__CUDACC_VER_MAJOR__ > 11 || ((__CUDACC_VER_MAJOR__ == 11) && (__CUDACC_VER_MINOR__ >= 2)))
#define AMREX_CUDA_GE_11_2 1
#endif

#if defined(AMREX_USE_HIP)
#define AMREX_HIP_OR_CUDA(a,b) a
#elif defined(AMREX_USE_CUDA)
#define AMREX_HIP_OR_CUDA(a,b) b
#else
#define AMREX_HIP_OR_CUDA(a,b) ((void)0);
#endif

#if defined(AMREX_USE_HIP)
#define AMREX_HIP_OR_CUDA_OR_SYCL(a,b,c) a
#elif defined(AMREX_USE_CUDA)
#define AMREX_HIP_OR_CUDA_OR_SYCL(a,b,c) b
#elif defined(AMREX_USE_SYCL)
#define AMREX_HIP_OR_CUDA_OR_SYCL(a,b,c) c
#else
#define AMREX_HIP_OR_CUDA_OR_SYCL(a,b,c) ((void)0);
#endif

#ifdef AMREX_USE_GPU
#define AMREX_GPU_OR_CPU(a,b) a
#else
#define AMREX_GPU_OR_CPU(a,b) b
#endif

#ifdef AMREX_USE_SYCL
#define AMREX_SYCL_ONLY(a) a
#else
#define AMREX_SYCL_ONLY(a) ((void)0)
#endif

#ifdef AMREX_USE_SYCL
#if (AMREX_SPACEDIM == 1)
#  define AMREX_SYCL_1D_ONLY(a) a
#  define AMREX_SYCL_2D_ONLY(a) ((void)0)
#  define AMREX_SYCL_3D_ONLY(a) ((void)0)
#elif (AMREX_SPACEDIM == 2)
#  define AMREX_SYCL_1D_ONLY(a) ((void)0)
#  define AMREX_SYCL_2D_ONLY(a) a
#  define AMREX_SYCL_3D_ONLY(a) ((void)0)
#elif (AMREX_SPACEDIM == 3)
#  define AMREX_SYCL_1D_ONLY(a) ((void)0)
#  define AMREX_SYCL_2D_ONLY(a) ((void)0)
#  define AMREX_SYCL_3D_ONLY(a) a
#endif
#else
#  define AMREX_SYCL_1D_ONLY(a) ((void)0)
#  define AMREX_SYCL_2D_ONLY(a) ((void)0)
#  define AMREX_SYCL_3D_ONLY(a) ((void)0)
#endif

namespace amrex {
    enum struct RunOn { Gpu, Cpu, Device=Gpu, Host=Cpu };
}

namespace amrex { // NOLINT(modernize-concat-nested-namespaces)

#ifdef AMREX_USE_HIP
using gpuStream_t = hipStream_t;
#elif defined(AMREX_USE_CUDA)
using gpuStream_t = cudaStream_t;
#endif

namespace Gpu {

#if defined(AMREX_USE_GPU)

    extern bool in_launch_region;

    [[nodiscard]] inline bool inLaunchRegion () noexcept { return in_launch_region; }
    [[nodiscard]] inline bool notInLaunchRegion () noexcept { return !in_launch_region; }

    /**
     * Enable/disable GPU kernel launches.
     *
     * \note This will only switch from GPU to CPU for kernels launched
     * with macros. Functions like `amrex::ParallelFor` will be unaffected.
     * Therefore it should not be used for comparing GPU to non-GPU performance
     * or behavior.
     *
     * \code
     *  Gpu::setLaunchRegion(0);
     *
     *  //...
     *
     *  Gpu::setLaunchRegion(1);
     * \endcode
     *
     * Will disable the launching of GPU kernels between the calls.
     */

    inline bool setLaunchRegion (bool launch) noexcept {
        bool r = in_launch_region;
        in_launch_region = launch;
        return r;
    }

    extern bool in_graph_region;
    [[nodiscard]] inline bool inGraphRegion() { return (in_graph_region && in_launch_region); }
    [[nodiscard]] inline bool notInGraphRegion() { return (!in_graph_region || !in_launch_region); }

    inline bool setGraphRegion (bool graph) {
        bool r = in_graph_region;
        in_graph_region = graph;
        return r;
    }

    struct [[nodiscard]] LaunchSafeGuard
    {
        explicit LaunchSafeGuard (bool flag) noexcept
            : m_old(setLaunchRegion(flag)) {}
        ~LaunchSafeGuard () { setLaunchRegion(m_old); }
    private:
        bool m_old;
    };

    struct [[nodiscard]] GraphSafeGuard
    {
        explicit GraphSafeGuard (bool flag) noexcept
            : m_old(setGraphRegion(flag)) {}
        ~GraphSafeGuard () { setGraphRegion(m_old); }
    private:
        bool m_old;
    };

    extern bool in_single_stream_region;
    extern bool in_nosync_region;

    [[nodiscard]] inline bool inSingleStreamRegion () noexcept { return in_single_stream_region; }
    [[nodiscard]] inline bool inNoSyncRegion () noexcept { return in_nosync_region; }

    inline bool setSingleStreamRegion (bool b) noexcept {
        return std::exchange(in_single_stream_region, b);
    }

    inline bool setNoSyncRegion (bool b) noexcept {
        return std::exchange(in_nosync_region, b);
    }

    /**
     * This struct provides a RAII-style mechansim for changing the number
     * of streams returned by Gpu::numStreams() to a single stream.
     */
    struct [[nodiscard]] SingleStreamRegion
    {
        SingleStreamRegion () noexcept
            : m_prev_flag(std::exchange(in_single_stream_region,true))
        {}

        ~SingleStreamRegion () { in_single_stream_region = m_prev_flag; }

    private:
        bool m_prev_flag;
    };

    /**
     * This struct provides a RAII-style mechansim for disabling GPU
     * synchronization in MFITer by default.  Note that explicit calls to
     * Gpu::steramSynchronize and Gpu::deviceSynchronize still work.
     */
    struct [[nodiscard]] NoSyncRegion
    {
        NoSyncRegion () noexcept
            : m_prev_flag(std::exchange(in_nosync_region,true))
        {}

        ~NoSyncRegion () { in_nosync_region = m_prev_flag; }

    private:
        bool m_prev_flag;
    };

#else

    [[nodiscard]] inline static constexpr bool inLaunchRegion () { return false; }
    [[nodiscard]] inline static constexpr bool notInLaunchRegion () { return true; }
    [[nodiscard]] inline static constexpr bool setLaunchRegion (bool) { return false; }

    [[nodiscard]] inline static constexpr bool inGraphRegion () { return false; }
    [[nodiscard]] inline static constexpr bool notInGraphRegion () { return true; }
    [[nodiscard]] inline static constexpr bool setGraphRegion (bool) { return false; }

    struct [[nodiscard]] LaunchSafeGuard
    {
        explicit LaunchSafeGuard (bool) {}
    };

    struct [[nodiscard]] GraphSafeGuard
    {
        explicit GraphSafeGuard (bool) {}
    };

    [[nodiscard]] inline static constexpr bool inSingleStreamRegion () { return false; }
    [[nodiscard]] inline static constexpr bool inNoSyncRegion () { return true; }
    [[nodiscard]] inline static constexpr bool setSingleStreamRegion (bool) { return false; }
    [[nodiscard]] inline static constexpr bool setNoSyncRegion (bool) { return true; }
    struct [[nodiscard]] SingleStreamRegion {};
    struct [[nodiscard]] NoSyncRegion {};

#endif

}
}

#endif
