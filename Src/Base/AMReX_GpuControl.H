#ifndef AMREX_GPU_CONTROL_H_
#define AMREX_GPU_CONTROL_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuTypes.H>

#ifndef AMREX_GPU_MAX_THREADS
#define AMREX_GPU_MAX_THREADS 256
#endif

#if defined(AMREX_USE_HIP)
#define AMREX_HIP_OR_CUDA(a,b) a
#elif defined(AMREX_USE_CUDA)
#define AMREX_HIP_OR_CUDA(a,b) b
#else
#define AMREX_HIP_OR_CUDA(a,b) ((void)0);
#endif

#if defined(AMREX_USE_HIP)
#define AMREX_HIP_OR_CUDA_OR_DPCPP(a,b,c) a
#elif defined(AMREX_USE_CUDA)
#define AMREX_HIP_OR_CUDA_OR_DPCPP(a,b,c) b
#elif defined(AMREX_USE_DPCPP)
#define AMREX_HIP_OR_CUDA_OR_DPCPP(a,b,c) c
#else
#define AMREX_HIP_OR_CUDA_OR_DPCPP(a,b,c) ((void)0);
#endif

#ifdef AMREX_USE_GPU
#define AMREX_GPU_OR_CPU(a,b) a
#else
#define AMREX_GPU_OR_CPU(a,b) b
#endif

#ifdef AMREX_USE_DPCPP
#define AMREX_DPCPP_ONLY(a) a
#else
#define AMREX_DPCPP_ONLY(a) ((void)0)
#endif

#ifdef AMREX_USE_DPCPP
#if (AMREX_SPACEDIM == 1)
#  define AMREX_DPCPP_1D_ONLY(a) a
#  define AMREX_DPCPP_2D_ONLY(a) ((void)0)
#  define AMREX_DPCPP_3D_ONLY(a) ((void)0)
#elif (AMREX_SPACEDIM == 2)
#  define AMREX_DPCPP_1D_ONLY(a) ((void)0)
#  define AMREX_DPCPP_2D_ONLY(a) a
#  define AMREX_DPCPP_3D_ONLY(a) ((void)0)
#elif (AMREX_SPACEDIM == 3)
#  define AMREX_DPCPP_1D_ONLY(a) ((void)0)
#  define AMREX_DPCPP_2D_ONLY(a) ((void)0)
#  define AMREX_DPCPP_3D_ONLY(a) a
#endif
#else
#  define AMREX_DPCPP_1D_ONLY(a) ((void)0)
#  define AMREX_DPCPP_2D_ONLY(a) ((void)0)
#  define AMREX_DPCPP_3D_ONLY(a) ((void)0)
#endif

namespace amrex {
    enum struct RunOn { Gpu, Cpu, Device=Gpu, Host=Cpu };
}

namespace amrex {

#ifdef AMREX_USE_HIP
using gpuStream_t = hipStream_t;
#elif defined(AMREX_USE_CUDA)
using gpuStream_t = cudaStream_t;
#endif

namespace Gpu {

#if defined(AMREX_USE_GPU)

    extern bool in_launch_region;

    inline bool inLaunchRegion() noexcept { return in_launch_region; }
    inline bool notInLaunchRegion() noexcept { return !in_launch_region; }

    inline bool setLaunchRegion (bool launch) noexcept {
        bool r = in_launch_region;
        in_launch_region = launch;
        return r;
    }

    extern bool in_graph_region;
    inline bool inGraphRegion() { return (in_graph_region && in_launch_region); }
    inline bool notInGraphRegion() { return (!in_graph_region || !in_launch_region); }

    inline bool setGraphRegion (bool graph) {
        bool r = in_graph_region;
        in_graph_region = graph;
        return r;
    }

    struct LaunchSafeGuard
    {
        explicit LaunchSafeGuard (bool flag) noexcept
            : m_old(setLaunchRegion(flag)) {}
        ~LaunchSafeGuard () { setLaunchRegion(m_old); }
    private:
        bool m_old;
    };

    struct GraphSafeGuard
    {
        explicit GraphSafeGuard (bool flag) noexcept
            : m_old(setGraphRegion(flag)) {}
        ~GraphSafeGuard () { setGraphRegion(m_old); }
    private:
        bool m_old;
    };

    struct ScopedDefaultStream
    {
    public:
        ScopedDefaultStream () noexcept;
        ~ScopedDefaultStream ();
    private:
        gpuStream_t m_prev_stream;
    };

#else

    inline static constexpr bool inLaunchRegion () { return false; }
    inline static constexpr bool notInLaunchRegion () { return true; }
    inline static constexpr bool setLaunchRegion (bool) { return false; }

    inline static constexpr bool inGraphRegion () { return false; }
    inline static constexpr bool notInGraphRegion () { return true; }
    inline static constexpr bool setGraphRegion (bool) { return false; }

    struct LaunchSafeGuard
    {
        explicit LaunchSafeGuard (bool) {}
    };

    struct GraphSafeGuard
    {
        explicit GraphSafeGuard (bool) {}
    };

    struct ScopedDefaultStream {};

#endif

}
}

#endif
