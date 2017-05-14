
#include <AMReX_BArena.H>

void*
amrex::BArena::alloc (std::size_t _sz)
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
amrex::BArena::free (void* pt)
{
#if (defined(CUDA) && defined(CUDA_UM))
    cudaFree(pt);
#else
    operator delete(pt);
#endif
}
