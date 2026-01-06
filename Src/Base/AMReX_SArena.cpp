#include <AMReX_SArena.H>
#include <AMReX_GpuDevice.H>

void*
amrex::SArena::alloc (std::size_t sz)
{
    void* pt = The_Arena()->alloc(sz);
    m_profiler.profile_alloc(pt, sz);
    return pt;
}

void
amrex::SArena::free (void* pt)
{
    m_profiler.profile_free(pt);
    amrex::Gpu::freeAsync(The_Arena(), pt);
}

bool
amrex::SArena::isDeviceAccessible () const
{
    return The_Arena()->isDeviceAccessible();
}

bool
amrex::SArena::isHostAccessible () const
{
    return The_Arena()->isHostAccessible();
}

bool
amrex::SArena::isManaged () const
{
    return The_Arena()->isManaged();
}

bool
amrex::SArena::isDevice () const
{
    return The_Arena()->isDevice();
}

bool
amrex::SArena::isPinned () const
{
    return The_Arena()->isPinned();
}
