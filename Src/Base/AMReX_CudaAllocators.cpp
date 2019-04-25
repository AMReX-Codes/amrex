#include <AMReX_CudaAllocators.H>

#ifdef AMREX_USE_CUDA
namespace amrex
{

namespace
{
    ThrustManagedAllocator<char> g_cached_allocator;
}

namespace Cuda
{
    ThrustManagedAllocator<char>& The_ThrustCachedAllocator () { return g_cached_allocator; };
}

}
#endif
