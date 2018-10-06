#ifndef AMREX_CUDA_UTILITY_H_
#define AMREX_CUDA_UTILITY_H_

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

#define AMREX_CUDA_HOST __host__
#define AMREX_CUDA_DEVICE __device__
#define AMREX_CUDA_GLOBAL __global__
#define AMREX_CUDA_HOST_DEVICE __host__ __device__

#else

#define AMREX_CUDA_HOST
#define AMREX_CUDA_DEVICE
#define AMREX_CUDA_HOST_DEVICE
#define AMREX_CUDA_GLOBAL

#endif


#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

// ************************************************
// CUDA versions

#define AMREX_SIMPLE_LAUNCH(numBlocks, numThreads, function, ...)  \
    { \
      function<<<numBlocks, numThreads, 0, amrex::Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define AMREX_SIMPLE_L_LAUNCH(runOn, numBlocks, numThreads, ...)  \
    { \
      if (runOn == RunOn::GPU) \
      { \
         amrex::launch_global<<<numBlocks, numThreads, 0, amrex::Device::cudaStream()>>>(__VA_ARGS__); \
         CudaErrorCheck(); \
      } \
      else \
      { \
         amrex::launch_host(__VA_ARGS__); \
      } \
    }

#define AMREX_BOX_LAUNCH(box, function, ...)  \
    { \
      dim3 numBlocks, numThreads; \
      amrex::Device::c_threads_and_blocks(box.loVect(), box.hiVect(), numBlocks, numThreads); \
      function<<<numBlocks, numThreads, 0, amrex::Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define AMREX_BOX_L_LAUNCH(runOn, box, ...)  \
    { \
      if (runOn == RunOn::GPU) \
      { \
         dim3 numBlocks, numThreads; \
         amrex::Device::c_threads_and_blocks(box.loVect(), box.hiVect(), numBlocks, numThreads); \
         amrex::launch_global<<<numBlocks, numThreads, 0, amrex::Device::cudaStream()>>>(__VA_ARGS__); \
         CudaErrorCheck(); \
      } \
      else \
      { \
         amrex::launch_host(__VA_ARGS__); \
      } \
    }

#define AMREX_BOXCOMPS_LAUNCH(box, comps, function, ...)  \
    { \
      dim3 numBlocks, numThreads; \
      amrex::Device::c_comps_threads_and_blocks(box.loVect(), box.hiVect(), comps, numBlocks, numThreads); \
      function<<<numBlocks, numThreads, 0, amrex::Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define AMREX_BOXCOMPS_L_LAUNCH(runOn, box, comps, ...)  \
    { \
      if (runOn = RunOn::GPU) \
      { \
         dim3 numBlocks, numThreads; \
         amrex::Device::c_comps_threads_and_blocks(box.loVect(), box.hiVect(), comps, numBlocks, numThreads); \
         amrex::launch_global<<<numBlocks, numThreads, 0, amrex::Device::cudaStream()>>>(__VA_ARGS__); \
         CudaErrorCheck(); \
      } \
      else \
      { \
         amrex::launch_host(__VA_ARGS__); \
      } \
    }

#define AMREX_PARTICLES_LAUNCH(num_particles, function, ...) \
    { \
      int numThreads, numBlocks; \
      amrex::Device::particle_threads_and_blocks(num_particles, numThreads, numBlocks); \
      function<<<numBlocks, numThreads, 0, amrex::Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define AMREX_PARTICLES_L_LAUNCH(runon, num_particles, ...) \
    { \
      if (runon = RunOn::GPU) \
      { \
         int numThreads, numBlocks; \
         amrex::Device::particle_threads_and_blocks(num_particles, numThreads, numBlocks); \
         amrex::launch_global<<<numBlocks, numThreads, 0, amrex::Device::cudaStream()>>>(__VA_ARGS__); \
         CudaErrorCheck(); \
      } \
      else \
      { \
         amrex::launch_host(__VA_ARGS__); \
      } \
    }

#else

// ************************************************
// CPU versions

#define AMREX_SIMPLE_LAUNCH(numBlocks, numThreads, function, ...) function(__VA_ARGS__);
#define AMREX_BOX_LAUNCH(box, function, ...) function (__VA_ARGS__);
#define AMREX_BOXCOMPS_LAUNCH(box, comps, function, ...) function (__VA_ARGS__);
#define AMREX_PARTICLES_LAUNCH(num_particles, function, ...) function (__VA_ARGS__);

#define AMREX_SIMPLE_L_LAUNCH(runOn, numBlocks, numThreads, ...) amrex::launch_host(__VA_ARGS__);
#define AMREX_BOX_L_LAUNCH(runOn, box, ...) amrex::launch_host(__VA_ARGS__); 
#define AMREX_BOXCOMPS_L_LAUNCH(runOn, box, comps, function, ...) amrex::launch_host (__VA_ARGS__);
#define AMREX_PARTICLES_L_LAUNCH(runOn, num_particles, ...) amrex::launch_host (__VA_ARGS__);

#endif

namespace amrex {

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

// CPU variation

template<class L>
void launch_host (L f0) { f0(); }

template<class L, class... Lambdas>
void launch_host (L f0, Lambdas... fs) { f0(); launch_host(fs...); }

// ************************************************

}

namespace amrex {

    class Box;
    class IntVect;

    AMREX_CUDA_HOST_DEVICE
    Box getThreadBox (const Box& bx);

    AMREX_CUDA_HOST_DEVICE
    Box getThreadBox (const Box& bx, const IntVect& typ);

    AMREX_CUDA_HOST_DEVICE
    void getThreadComponentBox (const Box& bx, Box& tbx, int comp);

    AMREX_CUDA_HOST_DEVICE
    void getThreadIndex (int &index, int &size, const int num_particles);

}

extern "C" {
    void* amrex_gpu_malloc (std::size_t size);
    void amrex_gpu_free (void* p);
}

#endif
