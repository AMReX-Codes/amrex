#ifndef AMREX_CUDA_MANAGED_H_
#define AMREX_CUDA_MANAGED_H_

#include <AMReX_CudaQualifiers.H>

namespace amrex {

struct Managed {

#ifdef AMREX_USE_CUDA

  void *operator new(size_t len)
  {
    void *ptr;
    cudaMallocManaged(&ptr, len);
    cudaDeviceSynchronize();
    return ptr;
  }

  void operator delete(void *ptr)
  {
    cudaDeviceSynchronize();
    cudaFree(ptr);
  }

#endif
};

struct Pinned {

#ifdef AMREX_USE_CUDA

  void *operator new(size_t len)
  {
    void *ptr;
    cudaMallocHost(&ptr, len);
    cudaDeviceSynchronize();
    return ptr;
  }

  void operator delete(void *ptr)
  {
    cudaDeviceSynchronize();
    cudaFreeHost(ptr);
  }

#endif
};

}  //namespace


#endif
