#ifndef AMREX_MANAGED_H_
#define AMREX_MANAGED_H_

namespace amrex {

// ************************************************
// CUDA versions

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

#define AMREX_CUDA_HOST __host__
#define AMREX_CUDA_DEVICE __device__
#define AMREX_CUDA_GLOBAL __global__

#define AMREX_SIMPLE_LAUNCH(numBlocks, numThreads, function, ...)  \
    { \
      function<<<numBlocks, numThreads, 0, Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define AMREX_SIMPLE_L_LAUNCH(numBlocks, numThreads, ...)  \
    { \
      launch_global<<<numBlocks, numThreads, 0, Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define AMREX_BOX_LAUNCH(box, function, ...)  \
    { \
      dim3 numBlocks, numThreads; \
      Device::c_threads_and_blocks(box.loVect(), box.hiVect(), numBlocks, numThreads); \
      function<<<numBlocks, numThreads, 0, Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define AMREX_BOXCOMPS_L_LAUNCH(box, comps, ...)  \
    { \
      dim3 numBlocks, numThreads; \
      Device::c_comps_threads_and_blocks(box.loVect(), box.hiVect(), comps, numBlocks, numThreads); \
      launch_global<<<numBlocks, numThreads, 0, Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define AMREX_BOXCOMPS_LAUNCH(box, comps, function, ...)  \
    { \
      dim3 numBlocks, numThreads; \
      Device::c_comps_threads_and_blocks(box.loVect(), box.hiVect(), comps, numBlocks, numThreads); \
      function<<<numBlocks, numThreads, 0, Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define AMREX_BOX_L_LAUNCH(box, ...)  \
    { \
      dim3 numBlocks, numThreads; \
      Device::c_threads_and_blocks(box.loVect(), box.hiVect(), numBlocks, numThreads); \
      launch_global<<<numBlocks, numThreads, 0, Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define AMREX_PARTICLES_LAUNCH(num_particles, function, ...) \
    { \
      int numThreads, numBlocks; \
      Device::particle_threads_and_blocks(num_particles, numThreads, numBlocks); \
      function<<<numBlocks, numThreads, 0, Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define AMREX_PARTICLES_L_LAUNCH(num_particles, ...) \
    { \
      int numThreads, numBlocks; \
      Device::particle_threads_and_blocks(num_particles, numThreads, numBlocks); \
      launch_global<<<numBlocks, numThreads, 0, Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

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
// CPU versions
#else

#define AMREX_SIMPLE_LAUNCH(numBlocks, numThreads, function, ...) function(__VA_ARGS__);
#define AMREX_BOX_LAUNCH(box, function, ...) function (__VA_ARGS__);
#define AMREX_PARTICLES_LAUNCH(num_particles, function, ...) function (__VA_ARGS__);

#define AMREX_SIMPLE_L_LAUNCH(numBlocks, numThreads, ...) launch_global(__VA_ARGS__);
#define AMREX_BOX_L_LAUNCH(box, ...) launch_global(__VA_ARGS__); 
#define AMREX_PARTICLES_L_LAUNCH(num_particles, ...) launch_global(__VA_ARGS__);

#define AMREX_CUDA_HOST
#define AMREX_CUDA_DEVICE
#define AMREX_CUDA_GLOBAL


inline bool syncDevice() {}
inline bool checkManaged(const void *ptr) {}
struct Managed {};
struct Pinned  {};

#endif

// ************************************************
//  Variadic lambda function wrappers for C++ CUDA Kernel calls.

template<class L>
AMREX_CUDA_GLOBAL void launch_global (L f0) { f0(); }

template<class L, class... Lambdas>
AMREX_CUDA_GLOBAL void launch_global (L f0, Lambdas... fs) { f0(); launch_device(fs...); }

template<class L>
AMREX_CUDA_DEVICE void launch_device (L f0) { f0(); }

template<class L, class... Lambdas>
AMREX_CUDA_DEVICE void launch_device (L f0, Lambdas... fs) { f0(); launch_device(fs...); }

// ************************************************

}  //namespace


#endif
