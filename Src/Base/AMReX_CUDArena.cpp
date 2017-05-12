
#include <AMReX_CUDArena.H>

#include <cuda_runtime.h>

void*
amrex::CUDArena::alloc (std::size_t _sz)
{
    void* pt;

    cudaMallocManaged((void**) &pt, _sz);

    return pt;
}

void
amrex::CUDArena::free (void* pt)
{
    cudaFree(pt);
}
