#include <AMReX.H>
#include <AMReX_Device.H>

#ifndef AMREX_MANAGED
#define AMREX_MANAGED

// ************************************************

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

#define AMREX_CUDA_HOST __host__
#define AMREX_CUDA_DEVICE __device__
#define AMREX_CUDA_GLOBAL __global__
/*
inline void CudaErrorCheck() {
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::string errStr("CUDA failure: ");
        errStr.append(cudaGetErrorString(err));
        amrex::Abort(errStr);
    }
}
*/
inline void syncDevice() 
{
   cudaDeviceSynchronize();
}

inline bool checkManaged(const void *ptr)
{
   cudaPointerAttributes ptr_attr;
   cudaPointerGetAttributes(&ptr_attr, ptr);
   return ptr_attr.isManaged;
}

#define AMREX_SIMPLE_LAUNCH(function, numBlocks, numThreads, ...)  \
    { \
      function<<<numBlocks, numThreads>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define AMREX_BOX_LAUNCH(box, function, ...)  \
    { \
      dim3 numBlocks, numThreads; \
      Device::c_threads_and_blocks(box.loVect(), box.hiVect(), numBlocks, numThreads); \
      function<<<numBlocks, numThreads>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define FORT_LAUNCH_PARTICLES(num_particles, function, ...) \
    { \
        int numThreads; int numBlocks; \
        Device::particle_threads_and_blocks(num_particles, numThreads, numBlocks); \
        function<<<numBlocks, numThreads>>>(__VA_ARGS__); \
        CudaErrorCheck(); \
    }



/*
Need particle_threads_and_blocks before turning this one on.

#define AMREX_PARTICLES_LAUNCH(num_particles, function, ...) \
    { \
      int numThreads, numBlocks; \
      Device::particle_threads_and_blocks(num_particles, numThreads, numBlocks); \
      function<<<numBlocks, numThreads, 0, Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }
*/
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

// ************************************************

#else

#define AMREX_SIMPLE_LAUNCH(function, numBlocks, numThreads, ...) function(__VA_ARGS__)
#define AMREX_BOX_LAUNCH(box, function, ...) function (__VA_ARGS__) 
#define AMREX_PARTICLES_LAUNCH(num_particles, function, ...) function (__VA_ARGS__)
#define FORT_LAUNCH_PARTICLES(num_particles, function, ...) function(__VA_ARGS__)
#define AMREX_CUDA_HOST
#define AMREX_CUDA_DEVICE
#define AMREX_CUDA_GLOBAL


inline bool syncDevice() {}
inline bool checkManaged(const void *ptr) {}
struct Managed { };
struct Pinned  { };

#endif

// ************************************************


#endif
