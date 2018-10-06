#ifndef AMREX_MANAGED_H_
#define AMREX_MANAGED_H_

#include <AMReX_CUDA_Utility.H>

namespace amrex {

#ifndef __CUDA_ARCH__

struct Managed {

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

};

struct Pinned {

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
};

#else

struct Managed { };
struct Pinned  { };

#endif


}  //namespace


#endif
