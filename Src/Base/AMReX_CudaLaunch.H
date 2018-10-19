#ifndef AMREX_CUDA_LAUNCH_H_
#define AMREX_CUDA_LAUNCH_H_

#include <AMReX_CudaQualifiers.H>
#include <cstddef>
#include <AMReX_Box.H>
#include <AMReX_CudaControl.H>
#include <AMReX_Device.H>

#define AMREX_CUDA_STRIDE 1
#define AMREX_CUDA_Y_STRIDE 1
#define AMREX_CUDA_Z_STRIDE 1

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

// ************************************************
// CUDA versions

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA(bbb,tbb,block) \
    if (amrex::Cuda::inLaunchRegion()) \
    { \
        auto amrex_i_st = Strategy(bbb); \
        amrex::launch_global<<<amrex_i_st.numBlocks, amrex_i_st.numThreads, 0, amrex::Device::cudaStream()>>>( \
        [=] AMREX_CUDA_DEVICE () { \
            long amrex_i_numpts = bbb.numPts(); \
            long amrex_i_inc = (amrex_i_numpts+AMREX_CUDA_STRIDE-1)/AMREX_CUDA_STRIDE; \
            for (long amrex_i_i = 0; amrex_i_i<amrex_i_numpts; amrex_i_i+=amrex_i_inc) \
            { \
                Box tbb = amrex::getThreadBox(bbb, amrex_i_i); \
                if (tbb.ok()) block \
            } \
        }); \
        CudaErrorCheck(); \
    } \
    else { \
        const amrex::Box& tbb = bbb; \
        block \
    }

#define AMREX_CUDA_LAUNCH_DEVICE_LAMBDA(bbb,tbb,block) \
    if (amrex::Cuda::inLaunchRegion()) \
    { \
        auto amrex_i_st = Strategy(bbb); \
        amrex::launch_global<<<amrex_i_st.numBlocks, amrex_i_st.numThreads, 0, amrex::Device::cudaStream()>>>( \
        [=] AMREX_CUDA_DEVICE () { \
            long amrex_i_numpts = bbb.numPts(); \
            long amrex_i_inc = (amrex_i_numpts+AMREX_CUDA_STRIDE-1)/AMREX_CUDA_STRIDE; \
            for (long amrex_i_i = 0; amrex_i_i<amrex_i_numpts; amrex_i_i+=amrex_i_inc) \
            { \
                Box tbb = amrex::getThreadBox(bbb, amrex_i_i); \
                if (tbb.ok()) block \
            } \
        }); \
        CudaErrorCheck(); \
    } \
    else { \
        amrex::Abort("AMREX_CUDA_LAUNCH_DEVICE_LAMBDA: cannot call device function from host"); \
    }

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_NOBOX(bbb,tbb,block) \
    if (amrex::Cuda::inLaunchRegion()) \
    { \
        auto amrex_i_st = Strategy(); \
        amrex::launch_global<<<amrex_i_st.numBlocks, amrex_i_st.numThreads, 0, amrex::Device::cudaStream()>>>( \
        [=] AMREX_CUDA_DEVICE () { \
            const auto amrex_i_lo = bbb.loVect3d(); \
            const auto amrex_i_hi = bbb.hiVect3d(); \
            for (int amrex_i_k = amrex_i_lo[2] + blockIdx.z * blockDim.z + threadIdx.z; amrex_i_k <= amrex_i_hi[2]; amrex_i_k += blockDim.z * gridDim.z) { \
            for (int amrex_i_j = amrex_i_lo[1] + blockIdx.y * blockDim.y + threadIdx.y; amrex_i_j <= amrex_i_hi[1]; amrex_i_j += blockDim.y * gridDim.y) { \
            for (int amrex_i_i = amrex_i_lo[0] + blockIdx.x * blockDim.x + threadIdx.x; amrex_i_i <= amrex_i_hi[0]; amrex_i_i += blockDim.x * gridDim.x) { \
                amrex::Box tbb(IntVect(AMREX_D_DECL(amrex_i_i,amrex_i_j,amrex_i_k)), \
                               IntVect(AMREX_D_DECL(amrex_i_i,amrex_i_j,amrex_i_k)), \
                               bbb.type()); \
                block \
            }}} \
        }); \
        CudaErrorCheck(); \
    } \
    else { \
        const amrex::Box& tbb = bbb; \
        block \
    }

#define AMREX_CUDA_LAUNCH_HOST_DEVICE(strategy, ...) \
    { \
      if (amrex::Cuda::inLaunchRegion()) \
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

#define AMREX_CUDA_LAUNCH_DEVICE(strategy, ...) \
    { \
      if (amrex::Cuda::inLaunchRegion()) \
      { \
         auto st = strategy; \
         amrex::launch_global<<<st.numBlocks, st.numThreads, 0, amrex::Device::cudaStream()>>>(__VA_ARGS__); \
         CudaErrorCheck(); \
      } \
      else \
      { \
         amrex::Abort("AMREX_CUDA_LAUNCH_DEVICE: cannot call device function from host"); \
      } \
    }

// Cannot respect Cuda::inLaunchRegion because function must be __global__.
#define AMREX_CUDA_LAUNCH_GLOBAL(strategy, function, ...) \
    { \
        auto st = strategy;                                             \
        function<<<st.numBlocks, st.numThreads, 0, amrex::Device::cudaStream()>>>(__VA_ARGS__); \
        CudaErrorCheck();                                               \
    }

#else

// ************************************************
// CPU versions

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA(bbb,tbb,block) \
    { \
        const amrex::Box& tbb = bbb; \
        block \
    }

#define AMREX_CUDA_LAUNCH_DEVICE_LAMBDA(bbb,tbb,block) \
    { \
        const amrex::Box& tbb = bbb; \
        block \
    }

#define AMREX_CUDA_LAUNCH_HOST_DEVICE(strategy, ...) amrex::launch_host(__VA_ARGS__);
#define AMREX_CUDA_LAUNCH_DEVICE(strategy, ...) amrex::launch_host(__VA_ARGS__);
#define AMREX_CUDA_LAUNCH_GLOBAL(strategy, function, ...) function(__VA_ARGS__);

#endif

namespace amrex {

// ************************************************
//  Variadic lambda function wrappers for C++ CUDA Kernel calls.

    template<class L>
    AMREX_CUDA_GLOBAL void launch_global (L f0) { f0(); }

    template<class L, class... Lambdas>
    AMREX_CUDA_GLOBAL void launch_global (L f0, Lambdas... fs) { f0(); call_device(fs...); }
    
    template<class L>
    AMREX_CUDA_DEVICE void call_device (L f0) { f0(); }
    
    template<class L, class... Lambdas>
    AMREX_CUDA_DEVICE void call_device (L f0, Lambdas... fs) { f0(); call_device(fs...); }
    
// CPU variation

    template<class L>
    void launch_host (L f0) { f0(); }
    
    template<class L, class... Lambdas>
    void launch_host (L f0, Lambdas... fs) { f0(); launch_host(fs...); }

// ************************************************

    struct ComponentBox {
        Box box;
        int ic;
        int nc;
    };

// ************************************************

    AMREX_CUDA_HOST_DEVICE
    Box getThreadBox (const Box& bx);

    AMREX_CUDA_HOST_DEVICE
    Box getThreadBox (const Box& bx, int n);

//    AMREX_CUDA_HOST_DEVICE
//    Box getThreadBox (const Box& bx, const IntVect& typ);

    AMREX_CUDA_HOST_DEVICE
    ComponentBox getThreadComponentBox (const Box& box, int ncomp);

    AMREX_CUDA_HOST_DEVICE
    void getThreadIndex (long &index, long &size, const long num_particles);

// ************************************************

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)
    struct Strategy {
        Strategy () {
            Device::grid_stride_threads_and_blocks(numBlocks,numThreads);
        }
        Strategy (const Box& box) {
            Device::n_threads_and_blocks( ((box.numPts()+AMREX_CUDA_STRIDE-1)/AMREX_CUDA_STRIDE), numBlocks, numThreads );
#if 0
            Box b = amrex::surroundingNodes(box);
            b -= box.smallEnd();
            b.coarsen(IntVect(AMREX_D_DECL(1,AMREX_CUDA_Y_STRIDE,AMREX_CUDA_Z_STRIDE)));
            Device::c_threads_and_blocks(b.loVect(), b.hiVect(), numBlocks, numThreads);
#endif
        }
        Strategy (const Box& box, int comps) {
            const Box& b = amrex::surroundingNodes(box);
            Device::c_comps_threads_and_blocks(b.loVect(), b.hiVect(), comps, numBlocks, numThreads);
        }
        Strategy (long N) {
            Device::n_threads_and_blocks(N, numBlocks, numThreads);
        }
        Strategy (dim3 nb, dim3 nt) : numBlocks(nb), numThreads(nt) {}
        
        dim3 numBlocks;
        dim3 numThreads;
    };
#endif

}

#endif
