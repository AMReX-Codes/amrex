#include <AMReX_BArena.H>
#if !defined(AMREX_FORTRAN_BOXLIB)
#include <AMReX_Gpu.H>
#endif

void*
amrex::BArena::alloc (std::size_t _sz)
{
    void* pt;

#if defined(AMREX_USE_GPU)
    if (device_use_managed_memory) {

	gpu_malloc_managed(&pt, &_sz);
	if (device_set_readonly)
	    Gpu::Device::mem_advise_set_readonly(pt, _sz);
	if (device_set_preferred) {
	    const int device = Gpu::Device::deviceId();
	    Gpu::Device::mem_advise_set_preferred(pt, _sz, device);
	}

    }
    else if (device_use_hostalloc) {

	gpu_hostalloc(&pt, &_sz);

    }
    else {

	gpu_malloc(&pt, &_sz);

    }
#else
    pt = ::operator new(_sz);
#endif

    return pt;
}

void
amrex::BArena::free (void* pt)
{
#if (defined(AMREX_USE_GPU)
    if (!device_use_hostalloc)
	gpu_free(pt);
    else
	gpu_freehost(pt);
#else
    ::operator delete(pt);
#endif
}

#ifdef AMREX_USE_GPU
void*
amrex::BArena::alloc_device (std::size_t _sz)
{
    void* pt = 0;

#ifdef AMREX_USE_GPU
    gpu_malloc(&pt, &_sz);
#endif

    return pt;
}

void
amrex::BArena::free_device (void* pt)
{
#ifdef AMREX_USE_GPU
    gpu_free(pt);
#endif
}
#endif
