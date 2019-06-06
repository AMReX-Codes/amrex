#ifndef AMREX_GPU_LAUNCH_H_
#define AMREX_GPU_LAUNCH_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuError.H>
#include <AMReX_GpuRange.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_Box.H>
#include <AMReX_Extension.H>
#include <AMReX_BLassert.H>
#include <AMReX_TypeTraits.H>
#include <cstddef>
#include <limits>

#ifdef AMREX_USE_GPU
#include <AMReX_GpuLaunchMacrosG.H>
#else
#include <AMReX_GpuLaunchMacrosC.H>
#endif

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_BOXIV(...) AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_BOXIV(__VA_ARGS__)
#define AMREX_LAUNCH_DEVICE_LAMBDA_BOXIV(...) AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_BOXIV(__VA_ARGS__)


#define AMREX_WRONG_NUM_ARGS(...) static_assert(false,"Wrong number of arguments to macro")
#define AMREX_GET_MACRO(_1,_2,_3,_4,_5,_6,_7,_8,_9,NAME,...) NAME
#define AMREX_LAUNCH_DEVICE_LAMBDA(...) AMREX_GET_MACRO(__VA_ARGS__,\
                                                        AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE_3, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE_2, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS)(__VA_ARGS__)

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA(...) AMREX_GET_MACRO(__VA_ARGS__,\
                                                        AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS)(__VA_ARGS__)


#define AMREX_FOR_1D(...)      AMREX_CUDA_DEVICE_FOR_1D(__VA_ARGS__)
#define AMREX_FOR_3D(...)      AMREX_CUDA_DEVICE_FOR_3D(__VA_ARGS__)
#define AMREX_FOR_4D(...)      AMREX_CUDA_DEVICE_FOR_4D(__VA_ARGS__)

#define AMREX_PARALLEL_FOR_1D(...)      AMREX_CUDA_DEVICE_PARALLEL_FOR_1D(__VA_ARGS__)
#define AMREX_PARALLEL_FOR_3D(...)      AMREX_CUDA_DEVICE_PARALLEL_FOR_3D(__VA_ARGS__)
#define AMREX_PARALLEL_FOR_4D(...)      AMREX_CUDA_DEVICE_PARALLEL_FOR_4D(__VA_ARGS__)

#define AMREX_HOST_DEVICE_FOR_1D(...) AMREX_CUDA_HOST_DEVICE_FOR_1D(__VA_ARGS__)
#define AMREX_HOST_DEVICE_FOR_3D(...) AMREX_CUDA_HOST_DEVICE_FOR_3D(__VA_ARGS__)
#define AMREX_HOST_DEVICE_FOR_4D(...) AMREX_CUDA_HOST_DEVICE_FOR_4D(__VA_ARGS__)


namespace amrex {

    // Because of CUDA, we cannot take rvalue lambdas.

// ************************************************
//  Variadic lambda function wrappers for C++ CUDA Kernel calls.

    template<class L>
    AMREX_GPU_GLOBAL void launch_global (L f0) { f0(); }

    template<class L, class... Lambdas>
    AMREX_GPU_GLOBAL void launch_global (L f0, Lambdas... fs) { f0(); call_device(fs...); }
    
    template<class L>
    AMREX_GPU_DEVICE void call_device (L f0) noexcept { f0(); }
    
    template<class L, class... Lambdas>
    AMREX_GPU_DEVICE void call_device (L f0, Lambdas... fs) noexcept { f0(); call_device(fs...); }
    
// CPU variation

    template<class L>
    void launch_host (L&& f0) noexcept { std::forward<L>(f0)(); }
    
    template<class L, class... Lambdas>
    void launch_host (L&& f0, Lambdas&&... fs) noexcept {
        std::forward<L>(f0)();
        launch_host(std::forward<Lambdas>(fs)...);
    }


    template <class T> class LayoutData;
    class FabArrayBase;

namespace Cuda {

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
    Box getThreadBox (const Box& bx, long offset) noexcept
    {
#if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
        const auto len = bx.length3d();
        long k = offset / (len[0]*len[1]);
        long j = (offset - k*(len[0]*len[1])) / len[0];
        long i = (offset - k*(len[0]*len[1])) - j*len[0];
        IntVect iv{AMREX_D_DECL(static_cast<int>(i),
                                static_cast<int>(j),
                                static_cast<int>(k))};
        iv += bx.smallEnd();
        return (bx & Box(iv,iv,bx.type()));
#else
        return bx;
#endif
    }

// ************************************************

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)
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
        ExecutionConfig (long N) noexcept {
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

#ifdef AMREX_USE_CUDA

template<typename T, typename L>
void launch (T const& n, L f) noexcept
{
    const auto ec = Cuda::ExecutionConfig(n);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        for (auto const i : Gpu::Range(n)) {
            f(i);
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void For (T n, L f) noexcept
{
    const auto ec = Cuda::ExecutionConfig(n);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        for (T i = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             i < n; i += stride) {
            f(i);
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void ParallelFor (T n, L f) noexcept
{
    For(n,f);
}

template <typename L>
void For (Box const& box, L f) noexcept
{
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    const auto ec = Cuda::ExecutionConfig(ncells);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, amrex::Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            f(i,j,k);
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename L>
void ParallelFor (Box const& box, L f) noexcept
{
    For(box,f);
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void For (Box const& box, T ncomp, L f) noexcept
{
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    const auto ec = Cuda::ExecutionConfig(ncells);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, amrex::Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            for (T n = 0; n < ncomp; ++n) {
                f(i,j,k,n);
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void ParallelFor (Box const& box, T ncomp, L f) noexcept
{
    For(box,ncomp,f);
}

#else

template<typename T, typename L>
AMREX_FORCE_INLINE
void launch (T const& n, L&& f) noexcept
{
    std::forward<L>(f)(n);
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
void For (T n, L&& f) noexcept
{
    for (T i = 0; i < n; ++i) {
        std::forward<L>(f)(i);
    }
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
void ParallelFor (T n, L&& f) noexcept
{
    AMREX_PRAGMA_SIMD
    for (T i = 0; i < n; ++i) {
        std::forward<L>(f)(i);
    }
}

template <typename L>
AMREX_FORCE_INLINE
void For (Box const& box, L&& f) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        std::forward<L>(f)(i,j,k);
    }}}
}

template <typename L>
AMREX_FORCE_INLINE
void ParallelFor (Box const& box, L&& f) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        std::forward<L>(f)(i,j,k);
    }}}
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
void For (Box const& box, T ncomp, L&& f) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (T n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            std::forward<L>(f)(i,j,k,n);
        }}}
    }
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
void ParallelFor (Box const& box, T ncomp, L&& f) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (T n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            std::forward<L>(f)(i,j,k,n);
        }}}
    }
}

#endif

}

#endif
