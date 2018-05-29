#include <iostream>
#include <AMReX_Print.H>
#include <AMReX_BArena.H>

void*
amrex::BArena::alloc (std::size_t _sz)
{
    amrex::Print() << "BArena::alloc(_sz)" << std::endl;
    return ::operator new(_sz);
}

void*
amrex::BArena::alloc (void* parent_ptr, std::size_t _sz)
{
  amrex::Print() << "BArena::alloc(ptr, _sz)" << std::endl;

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)  // Using CUDA

#ifdef __CUDA_ARCH__  // On Device
   void *ptr;
   cudaMalloc(&ptr, _sz);
   return ptr;

#else  // On Host
   cudaPointerAttributes ptr_attr;
   cudaPointerGetAttributes(&ptr_attr, parent_ptr);
   if (ptr_attr.isManaged)
   {
     void *ptr;
     cudaMallocManaged(&ptr, _sz);
     cudaDeviceSynchronize();
     return ptr;
   }
   else
   {
     return ::operator new(_sz);
   }
#endif

#else  // No CUDA
   return ::operator new(_sz);
#endif

}

void
amrex::BArena::free (void* ptr)
{
    amrex::Print() << "BArena::free()" << std::endl;

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)  // Using CUDA

#ifdef __CUDA_ARCH__  // On Device
   cudaFree(ptr);

#else  // On Host
   cudaPointerAttributes ptr_attr;
   cudaPointerGetAttributes(&ptr_attr, ptr);
   if (ptr_attr.isManaged)
   {
     cudaDeviceSynchronize();
     cudaFree(ptr); 
   }
   else
   {
     ::operator delete(ptr);
   }
#endif

#else  // No CUDA
   ::operator delete(ptr);
#endif
}
