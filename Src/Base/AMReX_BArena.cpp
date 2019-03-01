#include <AMReX_BArena.H>
#include <AMReX_Gpu.H>

void*
amrex::BArena::alloc (std::size_t sz_)
{
    return allocate_system(sz_);
}

void
amrex::BArena::free (void* pt)
{
    deallocate_system(pt);
}
