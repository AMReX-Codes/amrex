#ifndef AMREX_CUDA_LAUNCH_H_
#define AMREX_CUDA_LAUNCH_H_

#include <AMReX_CudaQualifiers.H>
#include <cstddef>
#include <AMReX_Box.H>
#include <AMReX_Device.H>

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

// ************************************************
// CUDA versions

#define AMREX_CUDA_LAUNCH_LAMBDA(strategy, ...) \
    { \
      if (amrex::Device::inLaunchRegion()) \
      { \
         auto st = strategy; \
         amrex::launch_global<<<st.numBlocks, st.numThreads, 0, amrex::Device::cudaStream()>>>(__VA_ARGS__); \
         CudaErrorCheck(); \
      } \
      else \
      { \
         amrex::launch_host(__VA_ARGS__); \
      } \
    }

// Cannot respect Device::inLaunchRegion because function must be __global__.
#define AMREX_CUDA_LAUNCH_GLOBAL(strategy, function, ...) \
    { \
        auto st = strategy;                                             \
        function<<<st.numBlocks, st.numThreads, 0, amrex::Device::cudaStream()>>>(__VA_ARGS__); \
        CudaErrorCheck();                                               \
    }

#else

// ************************************************
// CPU versions

#define AMREX_CUDA_LAUNCH_LAMBDA(strategy, ...) amrex::launch_host(__VA_ARGS__);
#define AMREX_CUDA_LAUNCH_GLOBAL(strategy, function, ...) function(__VA_ARGS__);

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

    AMREX_CUDA_HOST_DEVICE
    Box getThreadBox (const Box& bx);

//    AMREX_CUDA_HOST_DEVICE
//    Box getThreadBox (const Box& bx, const IntVect& typ);

    AMREX_CUDA_HOST_DEVICE
    void getThreadComponentBox (const Box& bx, Box& threadBox, int& scomp, int& ncomp);

    AMREX_CUDA_HOST_DEVICE
    void getThreadIndex (int &index, int &size, const int num_particles);

// ************************************************

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)
    struct Strategy {
//        Strategy () {
//            Device::grid_stride_threads_and_blocks(numBlocks,numThreads);
//        }
        Strategy (const Box& box) {
            const Box& b = amrex::surroundingNodes(box); // to ensure one thread per cell regardless box type
            Device::c_threads_and_blocks(b.loVect(), b.hiVect(), numBlocks, numThreads);
        }
        Strategy (const Box& box, int comps) {
            const Box& b = amrex::surroundingNodes(box);
            Device::c_comps_threads_and_blocks(b.loVect(), b.hiVect(), comps, numBlocks, numThreads);
        }
        Strategy (int N) {
            Device::n_threads_and_blocks(N, numBlocks, numThreads);
        }
        Strategy (dim3 nb, dim3 nt) : numBlocks(nb), numThreads(nt) {}
        
        dim3 numBlocks;
        dim3 numThreads;
    };
#endif

}

#endif
