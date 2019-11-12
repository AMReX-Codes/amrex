#include <AMReX_GpuAllocators.H>

#ifdef AMREX_USE_GPU
namespace amrex
{

namespace
{
    ThrustManagedAllocator<char> g_cached_allocator;
}

namespace Gpu 
{
    ThrustManagedAllocator<char>& The_ThrustCachedAllocator () { return g_cached_allocator; };
}

}
#endif
