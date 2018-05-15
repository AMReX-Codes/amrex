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

#define SIMPLE_LAUNCH(numBlocks, numThreads, function, ...)  \
    { \
      function<<<numBlocks, numThreads>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define DEVICE_FUNCTION __global__ 
#define DEVICE_CLASS_FUNCTION __host__ __device__
#define MANAGED public Managed

#else

struct Managed {

  void *operator new(size_t len)
  {
    void *ptr;
    return ptr;
  }

  void operator delete(void *ptr)
  {
    delete ptr;
  }

};

#define SIMPLE_LAUNCH(numBlocks, numThreads, function, ...) function(__VA_ARGS__)
#define DEVICE_FUNCTION
#define DEVICE_CLASS_FUNCTION 

#endif

// ************************************************


#endif
