#include <AMReX_BArena.H>
#ifdef AMREX_USE_DEVICE
#include <AMReX_Device.H>
#endif
#ifdef OMP_OFFLOAD
#include <cuda_runtime.h>
#endif

void*
amrex::BArena::alloc (std::size_t _sz)
{
    void* pt;

#if (defined(AMREX_USE_CUDA) && defined(CUDA_UM))
    if (device_use_managed_memory) {

	gpu_malloc_managed(&pt, &_sz);
	const int device = Device::deviceId();
	if (device_set_readonly)
	    Device::mem_advise_set_readonly(pt, _sz);
	if (device_set_preferred) {
	    const int device = Device::deviceId();
	    Device::mem_advise_set_preferred(pt, _sz, device);
	}

    }
    else if (device_use_hostalloc) {

	gpu_hostalloc(&pt, &_sz);

    }
    else {

	gpu_malloc(&pt, &_sz);

    }
// If combining the OpenMP offloading code with the existing CUDA_UM and
// AMREX_USE_CUDA code, the tutorials generate the wrong answer.  I do not know
// why. So until that gets sorted out, the OpenMP offloading code can't re-use
// any of the existing CUDA code here. This will lead to some duplication (such
// as these calls to cudaMallocManaged() and cudaFree()).
#elif (defined(OMP_OFFLOAD))
    cudaMallocManaged(&pt, _sz);
#else
    pt = operator new(_sz);
#endif

    return pt;
}

void
amrex::BArena::free (void* pt)
{
#if (defined(AMREX_USE_CUDA) && defined(CUDA_UM))
    if (!device_use_hostalloc)
	gpu_free(pt);
    else
	gpu_freehost(pt);
#elif (defined(OMP_OFFLOAD))
    cudaFree(pt);
#else
    operator delete(pt);
#endif
}

#ifdef AMREX_USE_DEVICE
void*
amrex::BArena::alloc_device (std::size_t _sz)
{
    void* pt = 0;

#ifdef AMREX_USE_CUDA
    gpu_malloc(&pt, &_sz);
#endif

    return pt;
}

void
amrex::BArena::free_device (void* pt)
{
#ifdef AMREX_USE_CUDA
    gpu_free(pt);
#endif
}
#endif
