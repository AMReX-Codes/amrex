#include <iostream>
#include <AMReX_Print.H>
#include <AMReX_BArena.H>
#include <AMReX_Device.H>

void*
amrex::BArena::alloc (std::size_t _sz)
{

// Using CUDA, CUDA is available, object created on the host.
// Manage everything for now. Improve later for other applications.
#if defined(AMREX_USE_CUDA) && defined(__CUDACC__) && !defined(__CUDA_ARCH__)

     void *ptr;
     cudaMallocManaged(&ptr, _sz);
     cudaDeviceSynchronize();
     CudaErrorCheck();
     return ptr;

#else  // No CUDA or on the device. 
   return ::operator new(_sz);
#endif
}

void
amrex::BArena::free (void* ptr)
{
#if defined(AMREX_USE_CUDA) && defined(__CUDACC__) && !defined(__CUDA_ARCH__)

   cudaDeviceSynchronize();
   cudaFree(ptr); 
   CudaErrorCheck();

#else 
   ::operator delete(ptr);
#endif
}
