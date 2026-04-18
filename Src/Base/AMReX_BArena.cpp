#include <AMReX_BArena.H>

void*
amrex::BArena::alloc (std::size_t sz_)
{
    sz_ = Arena::align(sz_);
#ifndef _WIN32
    void* pt = std::aligned_alloc(align_size, sz_);
#else
    void* pt = _aligned_malloc(sz_, align_size);
#endif
    m_profiler.profile_alloc(pt, sz_);
    return pt;
}

void
amrex::BArena::free (void* pt)
{
    m_profiler.profile_free(pt);
#ifndef _WIN32
    std::free(pt);
#else
    _aligned_free(pt);
#endif
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
