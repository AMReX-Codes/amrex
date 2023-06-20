#ifndef AMREX_GPU_LAUNCH_H_
#define AMREX_GPU_LAUNCH_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuKernelInfo.H>
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
#include <AMReX_RandomEngine.H>
#include <AMReX_Algorithm.H>
#include <cstddef>
#include <limits>
#include <algorithm>
#include <utility>

#define AMREX_GPU_NCELLS_PER_THREAD 3
#define AMREX_GPU_Y_STRIDE 1
#define AMREX_GPU_Z_STRIDE 1

#ifdef AMREX_USE_CUDA
#  define AMREX_LAUNCH_KERNEL(MT, blocks, threads, sharedMem, stream, ... ) \
        amrex::launch_global<MT><<<blocks, threads, sharedMem, stream>>>(__VA_ARGS__)
#elif defined(AMREX_USE_HIP)
#  define AMREX_LAUNCH_KERNEL(MT, blocks, threads, sharedMem, stream, ... ) \
        hipLaunchKernelGGL(launch_global<MT>, blocks, threads, sharedMem, stream, __VA_ARGS__)
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
            // If we change this, we must make sure it doesn't break FabArrayUtility Reduce*,
            // which assumes the decomposition is 1D.
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

    template <int MT>
    ExecutionConfig
    makeExecutionConfig (Long N) noexcept
    {
        ExecutionConfig ec(dim3{}, dim3{});
        ec.numBlocks.x = (std::max(N,Long(1)) + MT - 1) / MT;
        ec.numThreads.x = MT;
        AMREX_ASSERT(MT % Gpu::Device::warp_size == 0);
        return ec;
    }

    template <int MT>
    ExecutionConfig
    makeExecutionConfig (const Box& box) noexcept
    {
        ExecutionConfig ec(dim3{}, dim3{});
        ec.numBlocks.x = (std::max(box.numPts(),Long(1)) + MT - 1) / MT;
        ec.numThreads.x = MT;
        AMREX_ASSERT(MT % Gpu::Device::warp_size == 0);
        return ec;
    }
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

#include <AMReX_GpuLaunch.nolint.H>

#include <AMReX_CTOParallelForImpl.H>

#endif
