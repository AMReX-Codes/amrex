#include <AMReX_CudaAllocators.H>

namespace amrex
{

namespace
{
#ifdef AMREX_USE_CUDA
ThrustCachedAllocator g_allocator;
#endif
}

namespace Cuda
{
    ThrustCachedAllocator& The_ThrustCachedAllocator () { return g_allocator; };
}

}
