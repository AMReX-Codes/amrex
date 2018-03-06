#include <AMReX_BArena.H>

void*
amrex::BArena::alloc (std::size_t _sz)
{
    return ::operator new(_sz);
}

void
amrex::BArena::free (void* pt)
{
    ::operator delete(pt);
}
