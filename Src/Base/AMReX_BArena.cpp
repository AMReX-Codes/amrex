#include <AMReX_BArena.H>

void*
amrex::BArena::alloc (std::size_t sz_)
{
    return std::malloc(sz_);
}

void
amrex::BArena::free (void* pt)
{
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
