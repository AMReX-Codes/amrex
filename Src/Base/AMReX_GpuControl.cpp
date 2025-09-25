#include <AMReX_GpuControl.H>
#include <AMReX_BLassert.H>

#include <atomic>

namespace amrex::Gpu {

#if defined(AMREX_USE_GPU)
bool in_launch_region = true;
bool in_graph_region = false;
bool in_single_stream_region = false;
bool in_nosync_region = false;

namespace
{
    std::atomic<int> sync_launch_guard_counter{0};
}

void pushSyncLaunchGuard () noexcept
{
    sync_launch_guard_counter.fetch_add(1, std::memory_order_relaxed);
}

void popSyncLaunchGuard () noexcept
{
    int const previous = sync_launch_guard_counter.fetch_sub(1, std::memory_order_relaxed);
    AMREX_ASSERT(previous > 0);
    static_cast<void>(previous);
}

bool syncLaunchGuardActive () noexcept
{
    return sync_launch_guard_counter.load(std::memory_order_relaxed) > 0;
}

#else

void pushSyncLaunchGuard () noexcept {}
void popSyncLaunchGuard () noexcept {}
bool syncLaunchGuardActive () noexcept { return false; }

#endif

} // namespace amrex::Gpu
