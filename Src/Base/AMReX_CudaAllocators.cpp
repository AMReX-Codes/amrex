#include <AMReX_CudaAllocators.H>

#ifdef AMREX_USE_CUDA
namespace amrex
{

namespace
{
    ThrustDeviceAllocator<char> g_cached_allocator;
}

namespace Cuda
{
    ThrustDeviceAllocator<char>& The_ThrustCachedAllocator () { return g_cached_allocator; };
}

}
#endif
