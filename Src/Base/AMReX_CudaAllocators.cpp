#include <AMReX_CudaAllocators.H>

#ifdef AMREX_USE_CUDA
namespace amrex
{

namespace
{
    ThrustCachedAllocator g_allocator;
    ThrustPinnedAllocator g_pinned_allocator;
}

namespace Cuda
{
    ThrustCachedAllocator& The_ThrustCachedAllocator () { return g_allocator; };
    ThrustPinnedAllocator& The_ThrustPinnedAllocator () { return g_pinned_allocator; };
}

}
#endif
