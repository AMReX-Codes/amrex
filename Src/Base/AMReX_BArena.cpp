#include <AMReX_BArena.H>

void*
amrex::BArena::alloc (std::size_t sz_)
{
    void* pt = std::malloc(sz_);
    m_profiler.profile_alloc(pt, sz_);
    return pt;
}

void
amrex::BArena::free (void* pt)
{
    m_profiler.profile_free(pt);
    std::free(pt);
}

bool
amrex::BArena::isDeviceAccessible () const
{
    return false;
}

bool
amrex::BArena::isHostAccessible () const
{
    return true;
}

bool
amrex::BArena::isManaged () const
{
    return false;
}

bool
amrex::BArena::isDevice () const
{
    return false;
}

bool
amrex::BArena::isPinned () const
{
    return false;
}
