
#include <AMReX_CUDArena.H>

#include <cuda_runtime.h>

void*
amrex::CUDArena::alloc (std::size_t _sz)
{
    void* pt;

#if (defined(CUDA) && defined(CUDA_UM))
    cudaMallocManaged(&pt, _sz);
#else
    pt = operator new(_sz);
#endif

    return pt;
}

void
amrex::CUDArena::free (void* pt)
{
#if (defined(CUDA) && defined(CUDA_UM))
    cudaFree(pt);
#else
    operator delete(pt);
#endif
}
