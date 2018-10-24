#include <AMReX_CudaAllocators.H>

namespace amrex
{

namespace
{
#ifdef AMREX_USE_CUDA
ThrustCachedAllocator g_allocator;
#endif
}

#ifdef AMREX_USE_CUDA
namespace Cuda
{
    ThrustCachedAllocator& The_ThrustCachedAllocator () { return g_allocator; };
}
#endif

}
