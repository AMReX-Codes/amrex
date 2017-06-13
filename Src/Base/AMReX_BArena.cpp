
#include <AMReX_BArena.H>
#include <AMReX_Device.H>

void*
amrex::BArena::alloc (std::size_t _sz)
{
    void* pt;

#if (defined(CUDA) && defined(CUDA_UM))
    gpu_malloc_managed(&pt, &_sz);
    const int device = Device::cudaDeviceId();
    if (device_set_readonly)
	mem_advise_set_readonly(&pt, &_sz);
    if (device_set_preferred) {
	const int device = Device::cudaDeviceId();
	mem_advise_set_preferred(&pt, &_sz, &device);
    }
#else
    pt = operator new(_sz);
#endif

    return pt;
}

void
amrex::BArena::free (void* pt)
{
#if (defined(CUDA) && defined(CUDA_UM))
    gpu_free(pt);
#else
    operator delete(pt);
#endif
}

void*
amrex::BArena::alloc_device (std::size_t _sz)
{
    void* pt = 0;

#ifdef CUDA
    gpu_malloc(&pt, &_sz);
#endif

    return pt;
}

void
amrex::BArena::free_device (void* pt)
{
#ifdef CUDA
    gpu_free(pt);
#endif
}
