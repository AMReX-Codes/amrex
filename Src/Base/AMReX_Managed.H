#include <AMReX.H>

#ifndef AMREX_MANAGED
#define AMREX_MANAGED

// ************************************************

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

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

inline void CudaErrorCheck() {
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::string errStr("CUDA failure: ");
        errStr.append(cudaGetErrorString(err));
        amrex::Abort(errStr);
    }
}

#define AMREX_SIMPLE_LAUNCH(numBlocks, numThreads, function, ...)  \
    { \
      function<<<numBlocks, numThreads>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define AMREX_CUDA_HOST __host__
#define AMREX_CUDA_DEVICE __device__
#define AMREX_CUDA_GLOBAL __global__

// ************************************************

#else

struct Managed { };

#define AMREX_SIMPLE_LAUNCH(numBlocks, numThreads, function, ...) function(__VA_ARGS__)
#define AMREX_CUDA_HOST
#define AMREX_CUDA_DEVICE
#define AMREX_CUDA_GLOBAL

#endif

// ************************************************


#endif
