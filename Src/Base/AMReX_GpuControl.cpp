
#include <AMReX_GpuControl.H>
#include <AMReX_OpenMP.H>

namespace amrex::Gpu {

#if defined(AMREX_USE_GPU)

namespace {
    /* Per-thread execution-scope flags.
     *
     * These used to be four process-global bools, which the RAII guards save
     * and restore. That is unsound as soon as more than one host thread runs
     * AMReX, and not only for a caller that touches them directly: AMReX opens
     * the guards itself inside ordinary bound APIs -- Gpu::NoSyncRegion in
     * ParticleContainer::addParticles(), Gpu::LaunchSafeGuard in Gpu::AnyOf()
     * behind MultiFab::contains_nan()/contains_inf(). Two Python threads
     * calling add_particles() interleave as save(false) / save(true) /
     * restore(false) / restore(true), leaving in_nosync_region stuck true for
     * the rest of the process, after which AMReX silently stops synchronizing
     * streams. No caller discipline avoids that, which is why these became
     * per-thread rather than being documented away.
     *
     * File scope rather than class data members because MSVC forbids dllexport
     * on a thread_local -- the same reason Gpu::Device::getStreamIndex() is out
     * of line.
     */
    thread_local bool tl_in_launch_region        = true;
    thread_local bool tl_in_graph_region         = false;
    thread_local bool tl_in_single_stream_region = false;
    thread_local bool tl_in_nosync_region        = false;

#ifdef AMREX_USE_OMP
    /* What an OpenMP team observes.
     *
     * A team is created after the guard enclosing it, so its members do not
     * inherit the opening thread's thread_local. Opening a guard outside a
     * parallel region and having the team inside honour it is long-established
     * behaviour, so the setters broadcast here and a query from inside a team
     * reads this instead. Note the `if (Gpu::notInLaunchRegion())` clauses on
     * AMReX's own `omp parallel` directives are evaluated by the encountering
     * thread, which is outside the team, so those read the thread_local.
     *
     * Written only by a thread that is not itself in a parallel region, and
     * read only from inside one. Two host threads that each open a guard and
     * then each open a team still share these: an AMReX_OMP=ON build driven
     * from several host threads is out of scope, the same caveat the
     * per-thread RNG carries.
     */
    bool team_in_launch_region        = true;
    bool team_in_graph_region         = false;
    bool team_in_single_stream_region = false;
    bool team_in_nosync_region        = false;
#endif
}

bool inLaunchRegion () noexcept
{
#ifdef AMREX_USE_OMP
    if (OpenMP::in_parallel()) { return team_in_launch_region; }
#endif
    return tl_in_launch_region;
}

bool setLaunchRegion (bool launch) noexcept
{
    bool const r = tl_in_launch_region;
    tl_in_launch_region = launch;
#ifdef AMREX_USE_OMP
    if (!OpenMP::in_parallel()) { team_in_launch_region = launch; }
#endif
    return r;
}

bool inGraphRegionFlag () noexcept
{
#ifdef AMREX_USE_OMP
    if (OpenMP::in_parallel()) { return team_in_graph_region; }
#endif
    return tl_in_graph_region;
}

bool setGraphRegion (bool graph) noexcept
{
    bool const r = tl_in_graph_region;
    tl_in_graph_region = graph;
#ifdef AMREX_USE_OMP
    if (!OpenMP::in_parallel()) { team_in_graph_region = graph; }
#endif
    return r;
}

bool inSingleStreamRegion () noexcept
{
#ifdef AMREX_USE_OMP
    if (OpenMP::in_parallel()) { return team_in_single_stream_region; }
#endif
    return tl_in_single_stream_region;
}

bool setSingleStreamRegion (bool b) noexcept
{
    bool const r = tl_in_single_stream_region;
    tl_in_single_stream_region = b;
#ifdef AMREX_USE_OMP
    if (!OpenMP::in_parallel()) { team_in_single_stream_region = b; }
#endif
    return r;
}

bool inNoSyncRegion () noexcept
{
#ifdef AMREX_USE_OMP
    if (OpenMP::in_parallel()) { return team_in_nosync_region; }
#endif
    return tl_in_nosync_region;
}

bool setNoSyncRegion (bool b) noexcept
{
    bool const r = tl_in_nosync_region;
    tl_in_nosync_region = b;
#ifdef AMREX_USE_OMP
    if (!OpenMP::in_parallel()) { team_in_nosync_region = b; }
#endif
    return r;
}

#endif

}
