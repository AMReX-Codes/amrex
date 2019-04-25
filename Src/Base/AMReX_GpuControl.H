#ifndef AMREX_GPU_CONTROL_H_
#define AMREX_GPU_CONTROL_H_

#ifndef AMREX_CUDA_MAX_THREADS
#define AMREX_CUDA_MAX_THREADS 256
#endif

namespace amrex {
namespace Cuda {

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
        LaunchSafeGuard (bool flag) noexcept
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
        LaunchSafeGuard (bool) {}
    };

#endif

}
}

#endif
