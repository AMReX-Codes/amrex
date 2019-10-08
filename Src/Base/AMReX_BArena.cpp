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
