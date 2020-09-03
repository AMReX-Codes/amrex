#ifndef AMREX_GPU_LAUNCH_H_
#define AMREX_GPU_LAUNCH_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuTypes.H>
#include <AMReX_GpuError.H>
#include <AMReX_GpuRange.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_GpuReduce.H>
#include <AMReX_Tuple.H>
#include <AMReX_Box.H>
#include <AMReX_Loop.H>
#include <AMReX_Extension.H>
#include <AMReX_BLassert.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_GpuLaunchGlobal.H>
#include <cstddef>
#include <limits>
#include <algorithm>
#include <utility>

#define AMREX_GPU_NCELLS_PER_THREAD 3
#define AMREX_GPU_Y_STRIDE 1
#define AMREX_GPU_Z_STRIDE 1

#ifdef AMREX_USE_CUDA
#  define AMREX_LAUNCH_KERNEL(blocks, threads, sharedMem, stream, ... ) \
        amrex::launch_global<<<blocks, threads, sharedMem, stream>>>(__VA_ARGS__);
#elif defined(AMREX_USE_HIP) 
#  define AMREX_LAUNCH_KERNEL(blocks, threads, sharedMem, stream, ... ) \
        hipLaunchKernelGGL(launch_global, blocks, threads, sharedMem, stream, __VA_ARGS__);
#endif


namespace amrex {

// We cannot take rvalue lambdas.
// ************************************************
//  Variadic lambda function wrappers for C++ CUDA/HIP Kernel calls.

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    template<class L, class... Lambdas>
    AMREX_GPU_GLOBAL void launch_global (L f0, Lambdas... fs) { f0(); call_device(fs...); }

    template<class L>
    AMREX_GPU_DEVICE void call_device (L&& f0) noexcept { f0(); }
    
    template<class L, class... Lambdas>
    AMREX_GPU_DEVICE void call_device (L&& f0, Lambdas&&... fs) noexcept {
        f0();
        call_device(std::forward<Lambdas>(fs)...);
    }
#endif
    
// CPU variation

    template<class L>
    void launch_host (L&& f0) noexcept { f0(); }
    
    template<class L, class... Lambdas>
    void launch_host (L&& f0, Lambdas&&... fs) noexcept {
        f0();
        launch_host(std::forward<Lambdas>(fs)...);
    }


    template <class T> class LayoutData;
    class FabArrayBase;

namespace Gpu {

#ifdef AMREX_USE_GPU
    static constexpr std::size_t numThreadsPerBlockParallelFor () {
        return AMREX_GPU_MAX_THREADS;
    }
#else
    static constexpr std::size_t numThreadsPerBlockParallelFor () { return 0; }
#endif

// ************************************************

    struct ComponentBox {
        Box box;
        int ic;
        int nc;
    };

    struct GridSize {
        int numBlocks;
        int numThreads;
        int globalBlockId;
    };

// ************************************************

    void getGridSize (FabArrayBase const& fa, int ngrow, LayoutData<GridSize>& gs, int& ntotblocks);

    AMREX_GPU_HOST_DEVICE
    inline
    Box getThreadBox (const Box& bx, Long offset) noexcept
    {
#if AMREX_DEVICE_COMPILE
        const auto len = bx.length3d();
        Long k = offset / (len[0]*len[1]);
        Long j = (offset - k*(len[0]*len[1])) / len[0];
        Long i = (offset - k*(len[0]*len[1])) - j*len[0];
        IntVect iv{AMREX_D_DECL(static_cast<int>(i),
                                static_cast<int>(j),
                                static_cast<int>(k))};
        iv += bx.smallEnd();
        return (bx & Box(iv,iv,bx.type()));
#else
        amrex::ignore_unused(offset);
        return bx;
#endif
    }

// ************************************************

#ifdef AMREX_USE_GPU
    struct ExecutionConfig {
        ExecutionConfig () noexcept {
            Gpu::Device::grid_stride_threads_and_blocks(numBlocks,numThreads);
        }
        ExecutionConfig (const Box& box) noexcept {
            // If we change this, we must make sure it doesn't break say getGridSize, which
            // assumes the decomposition is 1D.  FabArrayUtility Reduce* as well.
            Gpu::Device::n_threads_and_blocks( ((box.numPts()+AMREX_GPU_NCELLS_PER_THREAD-1)/AMREX_GPU_NCELLS_PER_THREAD), numBlocks, numThreads );
#if 0
            Box b = amrex::surroundingNodes(box);
            b -= box.smallEnd();
            b.coarsen(IntVect(AMREX_D_DECL(1,AMREX_GPU_Y_STRIDE,AMREX_GPU_Z_STRIDE)));
            Gpu::Device::c_threads_and_blocks(b.loVect(), b.hiVect(), numBlocks, numThreads);
#endif
        }
        ExecutionConfig (const Box& box, int comps) noexcept {
            const Box& b = amrex::surroundingNodes(box);
            Gpu::Device::c_comps_threads_and_blocks(b.loVect(), b.hiVect(), comps, numBlocks, numThreads);
        }
        ExecutionConfig (Long N) noexcept {
            Gpu::Device::n_threads_and_blocks(N, numBlocks, numThreads);
        }
        ExecutionConfig (dim3 nb, dim3 nt, std::size_t sm=0) noexcept
            : numBlocks(nb), numThreads(nt), sharedMem(sm) {}
        
        dim3 numBlocks;
        dim3 numThreads;
        std::size_t sharedMem = 0;
    };
#endif

}
}


#ifdef AMREX_USE_GPU
#include <AMReX_GpuLaunchMacrosG.H>
#include <AMReX_GpuLaunchFunctsG.H>
#else
#include <AMReX_GpuLaunchMacrosC.H>
#include <AMReX_GpuLaunchFunctsC.H>
#endif

#define AMREX_WRONG_NUM_ARGS(...) static_assert(false,"Wrong number of arguments to macro")
#define AMREX_GET_MACRO(_1,_2,_3,_4,_5,_6,_7,_8,_9,NAME,...) NAME
#define AMREX_LAUNCH_DEVICE_LAMBDA(...) AMREX_GET_MACRO(__VA_ARGS__,\
                                                        AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_3, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS)(__VA_ARGS__)

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA(...) AMREX_GET_MACRO(__VA_ARGS__,\
                                                        AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS)(__VA_ARGS__)


#if (AMREX_SPACEDIM == 1)
#define AMREX_LAUNCH_DEVICE_LAMBDA_DIM(a1,a2,a3,b1,b2,b3,c1,c2,c3)      AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE     (a1,a2,a2)
#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM(a1,a2,a3,b1,b2,b3,c1,c2,c3) AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE(a1,a2,a3)
#elif (AMREX_SPACEDIM == 2)
#define AMREX_LAUNCH_DEVICE_LAMBDA_DIM(a1,a2,a3,b1,b2,b3,c1,c2,c3)      AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2     (a1,a2,a3,b1,b2,b3)
#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM(a1,a2,a3,b1,b2,b3,c1,c2,c3) AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2(a1,a2,a3,b1,b2,b3)
#elif (AMREX_SPACEDIM == 3)
#define AMREX_LAUNCH_DEVICE_LAMBDA_DIM(...)      AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_3     (__VA_ARGS__)
#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM(...) AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3(__VA_ARGS__)
#endif

#define AMREX_FOR_1D(...)      AMREX_GPU_DEVICE_FOR_1D(__VA_ARGS__)
#define AMREX_FOR_3D(...)      AMREX_GPU_DEVICE_FOR_3D(__VA_ARGS__)
#define AMREX_FOR_4D(...)      AMREX_GPU_DEVICE_FOR_4D(__VA_ARGS__)

#define AMREX_PARALLEL_FOR_1D(...)      AMREX_GPU_DEVICE_PARALLEL_FOR_1D(__VA_ARGS__)
#define AMREX_PARALLEL_FOR_3D(...)      AMREX_GPU_DEVICE_PARALLEL_FOR_3D(__VA_ARGS__)
#define AMREX_PARALLEL_FOR_4D(...)      AMREX_GPU_DEVICE_PARALLEL_FOR_4D(__VA_ARGS__)

#define AMREX_HOST_DEVICE_FOR_1D(...) AMREX_GPU_HOST_DEVICE_FOR_1D(__VA_ARGS__)
#define AMREX_HOST_DEVICE_FOR_3D(...) AMREX_GPU_HOST_DEVICE_FOR_3D(__VA_ARGS__)
#define AMREX_HOST_DEVICE_FOR_4D(...) AMREX_GPU_HOST_DEVICE_FOR_4D(__VA_ARGS__)

#define AMREX_HOST_DEVICE_PARALLEL_FOR_1D(...) AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_1D(__VA_ARGS__)
#define AMREX_HOST_DEVICE_PARALLEL_FOR_3D(...) AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_3D(__VA_ARGS__)
#define AMREX_HOST_DEVICE_PARALLEL_FOR_4D(...) AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_4D(__VA_ARGS__)

#ifdef AMREX_USE_GPU

#define AMREX_HOST_DEVICE_PARALLEL_FOR_1D_FLAG(where_to_run,n,i,block) \
    {  using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE (amrex_i_inttype i) noexcept \
            block \
        ); \
    } \
    else { \
        AMREX_PRAGMA_SIMD \
        for (amrex_i_inttype i = 0; i < n; ++i) { \
            block \
        } \
    }}

#define AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(where_to_run,box,i,j,k,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::LoopConcurrentOnCpu(box, [=] (int i, int j, int k) noexcept \
            block \
        ); \
    }

#define AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(where_to_run,box,nc,i,j,k,n,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(box, nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::LoopConcurrentOnCpu(box, nc, [=] (int i, int j, int k, int n) noexcept \
            block \
        ); \
    }

#define AMREX_HOST_DEVICE_FOR_1D_FLAG(where_to_run,n,i,block) \
    {  using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE (amrex_i_inttype i) noexcept \
            block \
        ); \
    } \
    else { \
        for (amrex_i_inttype i = 0; i < n; ++i) { \
            block \
        } \
    }}

#define AMREX_HOST_DEVICE_FOR_3D_FLAG(where_to_run,box,i,j,k,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::LoopOnCpu(box, [=] (int i, int j, int k) noexcept \
            block \
        ); \
    }

#define AMREX_HOST_DEVICE_FOR_4D_FLAG(where_to_run,box,nc,i,j,k,n,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(box, nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::LoopOnCpu(box, nc, [=] (int i, int j, int k, int n) noexcept \
            block \
        ); \
    }

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG(where_to_run,box,tbox,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        AMREX_LAUNCH_DEVICE_LAMBDA(box,tbox,block); \
    } else { \
        auto tbox = box; \
        block; \
    }

#else

#define AMREX_HOST_DEVICE_PARALLEL_FOR_1D_FLAG(where_to_run,n,i,block) \
    {  using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    amrex::ignore_unused(where_to_run); \
    AMREX_PRAGMA_SIMD \
    for (amrex_i_inttype i = 0; i < n; ++i) { \
        block \
    }}

#define AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(where_to_run,box,i,j,k,block) \
    amrex::ignore_unused(where_to_run); \
    amrex::LoopConcurrentOnCpu(box, [=] (int i, int j, int k) noexcept \
        block \
    );

#define AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(where_to_run,box,nc,i,j,k,n,block) \
    amrex::ignore_unused(where_to_run); \
    amrex::LoopConcurrentOnCpu(box, nc, [=] (int i, int j, int k, int n) noexcept \
        block \
    );

#define AMREX_HOST_DEVICE_FOR_1D_FLAG(where_to_run,n,i,block) \
    {  using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    amrex::ignore_unused(where_to_run); \
    for (amrex_i_inttype i = 0; i < n; ++i) { \
        block \
    }}

#define AMREX_HOST_DEVICE_FOR_3D_FLAG(where_to_run,box,i,j,k,block) \
    amrex::ignore_unused(where_to_run); \
    amrex::LoopOnCpu(box, [=] (int i, int j, int k) noexcept \
        block \
    );

#define AMREX_HOST_DEVICE_FOR_4D_FLAG(where_to_run,box,nc,i,j,k,n,block) \
    amrex::ignore_unused(where_to_run); \
    amrex::LoopOnCpu(box, nc, [=] (int i, int j, int k, int n) noexcept \
        block \
    );

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG(where_to_run,box,tbox,block) \
    amrex::ignore_unused(where_to_run); \
    { \
        auto tbox = box; \
        block; \
    }

#endif

#endif
