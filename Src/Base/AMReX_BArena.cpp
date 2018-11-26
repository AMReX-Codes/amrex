#include <AMReX_BArena.H>
#if !defined(AMREX_FORTRAN_BOXLIB)
#include <AMReX_Gpu.H>
#endif

void*
amrex::BArena::alloc (std::size_t sz_)
{
    void* pt;

#if defined(AMREX_USE_CUdA)
    if (device_use_managed_memory) {

	AMREX_GPU_SAFE_CALL(cudaMallocManaged(&pt, sz_));
	if (device_set_readonly)
	    Gpu::Device::mem_advise_set_readonly(pt, sz_);
	if (device_set_preferred) {
	    const int device = Gpu::Device::deviceId();
	    Gpu::Device::mem_advise_set_preferred(pt, sz_, device);
	}

    }
    else if (device_use_hostalloc) {

	AMREX_GPU_SAFE_CALL(cudaHostAlloc(&pt, sz_, cudaHostAllocMapped));

    }
    else {

	AMREX_GPU_SAFE_CALL(cudaMalloc(&pt, sz_));

    }
#else
    pt = ::operator new(sz_);
#endif

    return pt;
}

void
amrex::BArena::free (void* pt)
{
#if defined(AMREX_USE_GPU)
    if (!device_use_hostalloc) {
	AMREX_GPU_SAFE_CALL(cudaFree(pt));
    } else {
	AMREX_GPU_SAFE_CALL(cudaFreeHost(pt));
    }
#else
    ::operator delete(pt);
#endif
}
