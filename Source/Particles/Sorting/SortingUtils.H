#ifndef WARPX_PARTICLES_SORTING_SORTINGUTILS_H_
#define WARPX_PARTICLES_SORTING_SORTINGUTILS_H_

#include <AMReX_Gpu.H>
#include <AMReX_CudaContainers.H>

// TODO: Add documentation
void fillWithConsecutiveIntegers( amrex::Gpu::ManagedDeviceVector<long>& v ) {
#ifdef AMREX_USE_GPU
    // On GPU: Use thrust
    thrust::sequence( v.begin(), v.end() );
#else
    // On CPU: Use std library
    std::iota( v.begin(), v.end(), 0L );
#endif
}

#endif // WARPX_PARTICLES_SORTING_SORTINGUTILS_H_
