#ifndef AMREX_GPU_CONTROL_H_
#define AMREX_GPU_CONTROL_H_

#ifndef AMREX_GPU_MAX_THREADS
#define AMREX_GPU_MAX_THREADS 256
#endif

#if defined(AMREX_USE_HIP)
#define AMREX_HIP_OR_CUDA(a,b) a
#elif defined(AMREX_USE_CUDA)
#define AMREX_HIP_OR_CUDA(a,b) b
#else
#define AMREX_HIP_OR_CUDA(a,b) ((void)0)
#endif

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
#define AMREX_DEVICE_COMPILE 1
#endif

namespace amrex {
namespace Gpu {

#if defined(AMREX_USE_GPU)

    extern bool in_launch_region;

    inline bool inLaunchRegion () noexcept { return in_launch_region; }
    inline bool notInLaunchRegion () noexcept { return !in_launch_region; }

    inline bool setLaunchRegion (bool launch) noexcept {
        bool r =  in_launch_region;
        in_launch_region = launch;
        return r;
    }

    struct LaunchSafeGuard
    {
        explicit LaunchSafeGuard (bool flag) noexcept
            : m_old(setLaunchRegion(in_launch_region && flag)) {}
        ~LaunchSafeGuard () { setLaunchRegion(m_old); }
    private:
        bool m_old;
    };

#else

    inline static constexpr bool inLaunchRegion () { return false; }
    inline static constexpr bool notInLaunchRegion () { return true; }
    inline static constexpr bool setLaunchRegion (bool) { return false; }

    struct LaunchSafeGuard
    {
        explicit LaunchSafeGuard (bool) {}
    };

#endif

}
}

#endif
