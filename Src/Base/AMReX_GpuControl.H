#ifndef AMREX_GPU_CONTROL_H_
#define AMREX_GPU_CONTROL_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuTypes.H>

#include <utility>

#if defined(AMREX_USE_HIP) || defined(AMREX_USE_CUDA)
#define AMREX_GPU_STREAM_ALLOC_SUPPORT 1
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
    enum struct RunOn { Gpu, Cpu, Device=Gpu, Host=Cpu }; // NOLINT(readability-enum-initial-value)
}

#ifdef AMREX_USE_GPU
#define AMREX_DEFAULT_RUNON // need explicit RunOn when compiling for Gpu
#else
#define AMREX_DEFAULT_RUNON =amrex::RunOn::Host // by default run on Host when compiling for Cpu
#endif

namespace amrex { // NOLINT(modernize-concat-nested-namespaces)

#ifdef AMREX_USE_HIP
using gpuStream_t = hipStream_t;
#elif defined(AMREX_USE_CUDA)
using gpuStream_t = cudaStream_t;
#endif

namespace Gpu {

#if defined(AMREX_USE_GPU)

    /* The execution-scope flags (launch, graph, single-stream, nosync) are
     * per-thread: one host thread's RAII guard must not change what another
     * one does. An OpenMP team still observes a guard opened outside it,
     * because the setters broadcast to the team. See AMReX_GpuControl.cpp for
     * why this could not be left to caller discipline -- AMReX opens these
     * guards itself inside ordinary APIs.
     *
     * Out of line rather than inline over a thread_local, because MSVC forbids
     * dllexport on a thread_local; same treatment as
     * Gpu::Device::getStreamIndex().
     */
    [[nodiscard]] bool inLaunchRegion () noexcept;
    [[nodiscard]] inline bool notInLaunchRegion () noexcept { return !inLaunchRegion(); }

    /**
     * Enable/disable GPU kernel launches.
     *
     * \note This will only switch from GPU to CPU for kernels launched
     * with macros. Functions like `amrex::ParallelFor` will be unaffected.
     * Therefore it should not be used for comparing GPU to non-GPU performance
     * or behavior.
     *
     * \code{.cpp}
     *  Gpu::setLaunchRegion(0);
     *
     *  //...
     *
     *  Gpu::setLaunchRegion(1);
     * \endcode
     *
     * Will disable the launching of GPU kernels between the calls.
     */

    bool setLaunchRegion (bool launch) noexcept;

    //! The graph flag on its own; inGraphRegion() also requires the launch flag.
    [[nodiscard]] bool inGraphRegionFlag () noexcept;
    [[nodiscard]] inline bool inGraphRegion() { return (inGraphRegionFlag() && inLaunchRegion()); }
    [[nodiscard]] inline bool notInGraphRegion() { return (!inGraphRegionFlag() || !inLaunchRegion()); }

    bool setGraphRegion (bool graph) noexcept;

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

    [[nodiscard]] bool inSingleStreamRegion () noexcept;
    [[nodiscard]] bool inNoSyncRegion () noexcept;

    bool setSingleStreamRegion (bool b) noexcept;
    bool setNoSyncRegion (bool b) noexcept;

    /**
     * This struct provides a RAII-style mechanism for changing the number
     * of streams returned by Gpu::numStreams() to a single stream.
     */
    struct [[nodiscard]] SingleStreamRegion
    {
        SingleStreamRegion () noexcept
            : m_prev_flag(setSingleStreamRegion(true))
        {}

        ~SingleStreamRegion () { setSingleStreamRegion(m_prev_flag); }

    private:
        bool m_prev_flag;
    };

    /**
     * This struct provides a RAII-style mechanism for disabling GPU
     * synchronization in MFIter by default.  Note that explicit calls to
     * Gpu::streamSynchronize and Gpu::synchronize still work.
     */
    struct [[nodiscard]] NoSyncRegion
    {
        NoSyncRegion () noexcept
            : m_prev_flag(setNoSyncRegion(true))
        {}

        ~NoSyncRegion () { setNoSyncRegion(m_prev_flag); }

        [[nodiscard]] bool inNoSyncRegionPreviously () const { return m_prev_flag; }

    protected:
        bool m_prev_flag;
    };

#else

    [[nodiscard]] constexpr bool inLaunchRegion () { return false; }
    [[nodiscard]] constexpr bool notInLaunchRegion () { return true; }
    [[nodiscard]] constexpr bool setLaunchRegion (bool) { return false; }

    [[nodiscard]] constexpr bool inGraphRegion () { return false; }
    [[nodiscard]] constexpr bool notInGraphRegion () { return true; }
    [[nodiscard]] constexpr bool setGraphRegion (bool) { return false; }

    struct [[nodiscard]] LaunchSafeGuard
    {
        explicit LaunchSafeGuard (bool) {}
    };

    struct [[nodiscard]] GraphSafeGuard
    {
        explicit GraphSafeGuard (bool) {}
    };

    [[nodiscard]] constexpr bool inSingleStreamRegion () { return false; }
    [[nodiscard]] constexpr bool inNoSyncRegion () { return true; }
    [[nodiscard]] constexpr bool setSingleStreamRegion (bool) { return false; }
    [[nodiscard]] constexpr bool setNoSyncRegion (bool) { return true; }
    struct [[nodiscard]] SingleStreamRegion {};
    struct [[nodiscard]] NoSyncRegion {
        static constexpr bool inNoSyncRegionPreviously () { return true; }
    };

#endif

}
}

#endif
