#ifndef AMREX_GPU_CONTROL_H_
#define AMREX_GPU_CONTROL_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_GpuControl.H
 * \brief GPU launch-region and stream-region control flags and RAII guards.
 */

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuTypes.H>

#include <utility>

#if defined(AMREX_USE_CUDA) && (__CUDACC_VER_MAJOR__ > 11 || ((__CUDACC_VER_MAJOR__ == 11) && (__CUDACC_VER_MINOR__ >= 2)))
#define AMREX_CUDA_GE_11_2 1
#endif

#if defined(AMREX_USE_HIP) || defined(AMREX_CUDA_GE_11_2)
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
    /// \ingroup amrex_execution
    //! Tag to specify where a kernel or operation should run.
    enum struct RunOn { Gpu, Cpu, Device=Gpu, Host=Cpu }; // NOLINT(readability-enum-initial-value)
}

#ifdef AMREX_USE_GPU
#define AMREX_DEFAULT_RUNON // need explicit RunOn when compiling for Gpu
#else
#define AMREX_DEFAULT_RUNON =amrex::RunOn::Host // by default run on Host when compiling for Cpu
#endif

namespace amrex { // NOLINT(modernize-concat-nested-namespaces)

/// \cond DOXYGEN_IGNORE
#ifdef AMREX_USE_HIP
using gpuStream_t = hipStream_t;
#elif defined(AMREX_USE_CUDA)
using gpuStream_t = cudaStream_t;
#endif
/// \endcond

namespace Gpu {

#if defined(AMREX_USE_GPU)

    extern bool in_launch_region;

    /// \ingroup amrex_execution
    //! Return true if GPU kernel launches are currently enabled.
    [[nodiscard]] inline bool inLaunchRegion () noexcept { return in_launch_region; }
    /// \ingroup amrex_execution
    //! Return true if GPU kernel launches are currently disabled.
    [[nodiscard]] inline bool notInLaunchRegion () noexcept { return !in_launch_region; }

    /**
     * \ingroup amrex_execution
     * \brief Enable or disable GPU kernel launches; returns the previous flag value.
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
    inline bool setLaunchRegion (bool launch) noexcept {
        bool r = in_launch_region;
        in_launch_region = launch;
        return r;
    }

    extern bool in_graph_region;
    /// \ingroup amrex_execution
    //! Return true if both the graph region and the launch region are active.
    [[nodiscard]] inline bool inGraphRegion() { return (in_graph_region && in_launch_region); }
    /// \ingroup amrex_execution
    //! Return true if the graph region or the launch region is inactive.
    [[nodiscard]] inline bool notInGraphRegion() { return (!in_graph_region || !in_launch_region); }

    /// \ingroup amrex_execution
    //! Enable or disable the CUDA graph capture region; returns the previous flag value.
    inline bool setGraphRegion (bool graph) {
        bool r = in_graph_region;
        in_graph_region = graph;
        return r;
    }

    /**
     * \class amrex::Gpu::LaunchSafeGuard
     * \ingroup amrex_execution
     * \brief RAII guard that sets the GPU launch region flag on construction
     *        and restores it on destruction.
     */
    /// RAII guard that toggles the launch-region flag.
    struct [[nodiscard]] LaunchSafeGuard
    {
        //! Construct; sets the launch region flag to \p flag.
        explicit LaunchSafeGuard (bool flag) noexcept
            : m_old(setLaunchRegion(flag)) {}
        //! Destructor; restores the original launch region flag.
        ~LaunchSafeGuard () { setLaunchRegion(m_old); }
    private:
        bool m_old;
    };

    /**
     * \class amrex::Gpu::GraphSafeGuard
     * \ingroup amrex_execution
     * \brief RAII guard that sets the GPU graph region flag on construction
     *        and restores it on destruction.
     */
    /// RAII guard that toggles the CUDA-graph capture flag.
    struct [[nodiscard]] GraphSafeGuard
    {
        //! Construct; sets the graph region flag to \p flag.
        explicit GraphSafeGuard (bool flag) noexcept
            : m_old(setGraphRegion(flag)) {}
        //! Destructor; restores the original graph region flag.
        ~GraphSafeGuard () { setGraphRegion(m_old); }
    private:
        bool m_old;
    };

    extern bool in_single_stream_region;
    extern bool in_nosync_region;

    /// \ingroup amrex_execution
    //! Return true if the single-stream region is active.
    [[nodiscard]] inline bool inSingleStreamRegion () noexcept { return in_single_stream_region; }
    /// \ingroup amrex_execution
    //! Return true if the no-sync region is active.
    [[nodiscard]] inline bool inNoSyncRegion () noexcept { return in_nosync_region; }

    /// \ingroup amrex_execution
    //! Enable or disable the single-stream region; returns the previous flag value.
    inline bool setSingleStreamRegion (bool b) noexcept {
        return std::exchange(in_single_stream_region, b);
    }

    /// \ingroup amrex_execution
    //! Enable or disable the no-sync region; returns the previous flag value.
    inline bool setNoSyncRegion (bool b) noexcept {
        return std::exchange(in_nosync_region, b);
    }

    /**
     * \class amrex::Gpu::SingleStreamRegion
     * \ingroup amrex_execution
     * \brief RAII guard that restricts Gpu::numStreams() to a single stream
     *        for the duration of its lifetime.
     */
    struct [[nodiscard]] SingleStreamRegion
    {
        //! Construct; activates the single-stream region.
        SingleStreamRegion () noexcept
            : m_prev_flag(std::exchange(in_single_stream_region,true))
        {}
        //! Destructor; restores the previous stream-count state.
        ~SingleStreamRegion () { in_single_stream_region = m_prev_flag; }

    private:
        bool m_prev_flag;
    };

    /**
     * \class amrex::Gpu::NoSyncRegion
     * \ingroup amrex_execution
     * \brief RAII guard that disables automatic GPU synchronization in MFIter
     *        for the duration of its lifetime.
     *
     * Explicit calls to Gpu::streamSynchronize and Gpu::deviceSynchronize
     * are still honoured inside a NoSyncRegion.
     */
    struct [[nodiscard]] NoSyncRegion
    {
        //! Construct; activates the no-sync region.
        NoSyncRegion () noexcept
            : m_prev_flag(std::exchange(in_nosync_region,true))
        {}
        //! Destructor; restores the previous sync state.
        ~NoSyncRegion () { in_nosync_region = m_prev_flag; }

    private:
        bool m_prev_flag;
    };

/// \cond DOXYGEN_IGNORE
#else

    [[nodiscard]] constexpr bool inLaunchRegion () { return false; }
    [[nodiscard]] constexpr bool notInLaunchRegion () { return true; }
    [[nodiscard]] constexpr bool setLaunchRegion (bool) { return false; }

    [[nodiscard]] constexpr bool inGraphRegion () { return false; }
    [[nodiscard]] constexpr bool notInGraphRegion () { return true; }
    [[nodiscard]] constexpr bool setGraphRegion (bool) { return false; }

    /// RAII guard that toggles the launch-region flag (CPU build stub).
    struct [[nodiscard]] LaunchSafeGuard
    {
        explicit LaunchSafeGuard (bool) {}
    };

    /// RAII guard that toggles the CUDA-graph capture flag (CPU build stub).
    struct [[nodiscard]] GraphSafeGuard
    {
        explicit GraphSafeGuard (bool) {}
    };

    [[nodiscard]] constexpr bool inSingleStreamRegion () { return false; }
    [[nodiscard]] constexpr bool inNoSyncRegion () { return true; }
    [[nodiscard]] constexpr bool setSingleStreamRegion (bool) { return false; }
    [[nodiscard]] constexpr bool setNoSyncRegion (bool) { return true; }
    struct [[nodiscard]] SingleStreamRegion {};
    struct [[nodiscard]] NoSyncRegion {};

#endif
/// \endcond

}
}

#endif
