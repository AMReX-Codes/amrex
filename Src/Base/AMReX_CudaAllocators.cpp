#include <AMReX_CudaAllocators.H>

#ifdef AMREX_USE_CUDA
namespace amrex
{

namespace
{
    ThrustCachedAllocator g_allocator;
}

namespace Cuda
{
    ThrustCachedAllocator& The_ThrustCachedAllocator () { return g_allocator; };
}

}
#endif
