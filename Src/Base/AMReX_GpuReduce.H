#ifndef AMREX_GPU_REDUCE_H_
#define AMREX_GPU_REDUCE_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuTypes.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuUtility.H>
#include <AMReX_Functional.H>

#if !defined(AMREX_USE_CUB) && defined(AMREX_USE_CUDA) && defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 11)
#define AMREX_USE_CUB 1
#endif

#ifdef AMREX_USE_CUB
#include <cub/cub.cuh>
#endif

//
// Public interface
//
namespace amrex { namespace Gpu {
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void deviceReduceSum (T * dest, T source, Gpu::Handler const& h) noexcept;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void deviceReduceMin (T * dest, T source, Gpu::Handler const& h) noexcept;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void deviceReduceMax (T * dest, T source, Gpu::Handler const& h) noexcept;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void deviceReduceLogicalAnd (int * dest, int source, Gpu::Handler const& h) noexcept;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void deviceReduceLogicalOr (int * dest, int source, Gpu::Handler const& h) noexcept;
}}

//
// Reduce functions based on _shfl_down_sync
//

namespace amrex { namespace Gpu {

#ifdef AMREX_USE_DPCPP

template <int warpSize, typename T, typename F>
struct warpReduce
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T operator() (T x, sycl::ONEAPI::sub_group const& sg) const noexcept
    {
        for (int offset = warpSize/2; offset > 0; offset /= 2) {
            T y = sg.shuffle_down(x, offset);
            x = F()(x,y);
        }
        return x;
    }
};

template <int warpSize, typename T, typename WARPREDUCE>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduce (T x, WARPREDUCE && warp_reduce, T x0, Gpu::Handler const& h)
{
    T* shared = (T*)h.local;
    int tid = h.item->get_local_id(0);
    sycl::ONEAPI::sub_group const& sg = h.item->get_sub_group();
    int lane = sg.get_local_id()[0];
    int wid = sg.get_group_id()[0];
    int numwarps = sg.get_group_range()[0];
    x = warp_reduce(x, sg);
    // __syncthreads() prior to writing to shared memory is necessary
    // if this reduction call is occurring multiple times in a kernel,
    // and since we don't know how many times the user is calling it,
    // we do it always to be safe.
    h.item->barrier(sycl::access::fence_space::local_space);
    if (lane == 0) { shared[wid] = x; }
    h.item->barrier(sycl::access::fence_space::local_space);
    bool b = (tid == 0) || (tid < numwarps);
    x =  b ? shared[lane] : x0;
    if (wid == 0) { x = warp_reduce(x, sg); }
    return x;
}

template <int warpSize, typename T, typename WARPREDUCE, typename ATOMICOP>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void blockReduce_partial (T* dest, T x, WARPREDUCE && warp_reduce, ATOMICOP && atomic_op,
                          Gpu::Handler const& handler)
{
   sycl::ONEAPI::sub_group const& sg = handler.item->get_sub_group();
   int wid = sg.get_group_id()[0];
   if ((wid+1)*warpSize <= handler.numActiveThreads) {
       x = warp_reduce(x, sg); // full warp
       if (sg.get_local_id()[0] == 0) { atomic_op(dest, x); }
   } else {
       atomic_op(dest, x);
   }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceSum_full (T * dest, T source, Gpu::Handler const& h) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Plus<T> >(), (T)0, h);
    if (h.item->get_local_id(0) == 0) { Gpu::Atomic::AddNoRet(dest, source); }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMin_full (T * dest, T source, Gpu::Handler const& h) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Less<T> >(), source, h);
    if (h.item->get_local_id(0) == 0) { Gpu::Atomic::Min(dest, source); }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMax_full (T * dest, T source, Gpu::Handler const& h) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Greater<T> >(), source, h);
    if (h.item->get_local_id(0) == 0) { Gpu::Atomic::Max(dest, source); }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalAnd_full (int * dest, int source, Gpu::Handler const& h) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalAnd<int> >(), 1, h);
    if (h.item->get_local_id(0) == 0) { Gpu::Atomic::LogicalAnd(dest, source); }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalOr_full (int * dest, int source, Gpu::Handler const& h) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalOr<int> >(), 0, h);
    if (h.item->get_local_id(0) == 0) { Gpu::Atomic::LogicalOr(dest, source); }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceSum (T * dest, T source, Gpu::Handler const& h) noexcept
{
    if (h.isFullBlock()) {
        deviceReduceSum_full(dest, source, h);
    } else {
        Gpu::blockReduce_partial<Gpu::Device::warp_size>
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Plus<T> >(),
             Gpu::AtomicAdd<T>(), h);
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMin (T * dest, T source, Gpu::Handler const& h) noexcept
{
    if (h.isFullBlock()) {
        deviceReduceMin_full(dest, source, h);
    } else {
        Gpu::blockReduce_partial<Gpu::Device::warp_size>
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Less<T> >(),
             Gpu::AtomicMin<T>(), h);
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMax (T * dest, T source, Gpu::Handler const& h) noexcept
{
    if (h.isFullBlock()) {
        deviceReduceMax_full(dest, source, h);
    } else {
        Gpu::blockReduce_partial<Gpu::Device::warp_size>
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Greater<T> >(),
             Gpu::AtomicMax<T>(), h);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalAnd (int * dest, int source, Gpu::Handler const& h) noexcept
{
    if (h.isFullBlock()) {
        deviceReduceLogicalAnd_full(dest, source, h);
    } else {
        Gpu::blockReduce_partial<Gpu::Device::warp_size>
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalAnd<int> >(),
             Gpu::AtomicLogicalAnd<int>(), h);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalOr (int * dest, int source, Gpu::Handler const& h) noexcept
{
    if (h.isFullBlock()) {
        deviceReduceLogicalOr_full(dest, source, h);
    } else {
        Gpu::blockReduce_partial<Gpu::Device::warp_size>
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalOr<int> >(),
             Gpu::AtomicLogicalOr<int>(), h);
    }
}

#elif defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

template <int warpSize, typename T, typename F>
struct warpReduce
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T operator() (T x) const noexcept
    {
        for (int offset = warpSize/2; offset > 0; offset /= 2) {
            AMREX_HIP_OR_CUDA(T y = __shfl_down(x, offset);,
                              T y = __shfl_down_sync(0xffffffff, x, offset); )
            x = F()(x,y);
        }
        return x;
    }
};

template <int warpSize, typename T, typename WARPREDUCE>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduce (T x, WARPREDUCE && warp_reduce, T x0)
{
    __shared__ T shared[warpSize];
    int lane = threadIdx.x % warpSize;
    int wid = threadIdx.x / warpSize;
    x = warp_reduce(x);
    // __syncthreads() prior to writing to shared memory is necessary
    // if this reduction call is occurring multiple times in a kernel,
    // and since we don't know how many times the user is calling it,
    // we do it always to be safe.
    __syncthreads();
    if (lane == 0) { shared[wid] = x; }
    __syncthreads();
    bool b = (threadIdx.x == 0) || (threadIdx.x < blockDim.x / warpSize);
    x =  b ? shared[lane] : x0;
    if (wid == 0) { x = warp_reduce(x); }
    return x;
}

template <int warpSize, typename T, typename WARPREDUCE, typename ATOMICOP>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void blockReduce_partial (T* dest, T x, WARPREDUCE && warp_reduce, ATOMICOP && atomic_op,
                          Gpu::Handler const& handler)
{
    int warp = (int)threadIdx.x / warpSize;
    if ((warp+1)*warpSize <= handler.numActiveThreads) {
        x = warp_reduce(x); // full warp
        if (threadIdx.x % warpSize == 0) { atomic_op(dest, x); }
    } else {
        atomic_op(dest,x);
    }
}

#if defined(AMREX_USE_CUB)
template <int BLOCKDIMX, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceSum_full (T * dest, T source) noexcept
{
    typedef cub::BlockReduce<T,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::TempStorage temp_storage;
    // __syncthreads() prior to writing to shared memory is necessary
    // if this reduction call is occurring multiple times in a kernel,
    // and since we don't know how many times the user is calling it,
    // we do it always to be safe.
    __syncthreads();
    // Compute the block-wide reduction for thread0
    T aggregate = BlockReduce(temp_storage).Sum(source);
    if (threadIdx.x == 0) { Gpu::Atomic::AddNoRet(dest, aggregate); }
}
#endif

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceSum_full (T * dest, T source) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Plus<T> >(), (T)0);
    if (threadIdx.x == 0) { Gpu::Atomic::AddNoRet(dest, source); }
}

#if defined(AMREX_USE_CUB)
template <int BLOCKDIMX, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMin_full (T * dest, T source) noexcept
{
    typedef cub::BlockReduce<T,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::TempStorage temp_storage;
    // __syncthreads() prior to writing to shared memory is necessary
    // if this reduction call is occurring multiple times in a kernel,
    // and since we don't know how many times the user is calling it,
    // we do it always to be safe.
    __syncthreads();
    // Compute the block-wide reduction for thread0
    T aggregate = BlockReduce(temp_storage).Reduce(source, cub::Min());
    if (threadIdx.x == 0) { Gpu::Atomic::Min(dest, aggregate); }
}
#endif

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMin_full (T * dest, T source) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Less<T> >(), source);
    if (threadIdx.x == 0) { Gpu::Atomic::Min(dest, source); }
}

#if defined(AMREX_USE_CUB)
template <int BLOCKDIMX, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMax_full (T * dest, T source) noexcept
{
    typedef cub::BlockReduce<T,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::TempStorage temp_storage;
    // __syncthreads() prior to writing to shared memory is necessary
    // if this reduction call is occurring multiple times in a kernel,
    // and since we don't know how many times the user is calling it,
    // we do it always to be safe.
    __syncthreads();
    // Compute the block-wide reduction for thread0
    T aggregate = BlockReduce(temp_storage).Reduce(source, cub::Max());
    if (threadIdx.x == 0) { Gpu::Atomic::Max(dest, aggregate); }
}
#endif

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMax_full (T * dest, T source) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Greater<T> >(), source);
    if (threadIdx.x == 0) { Gpu::Atomic::Max(dest, source); }
}

#if defined(AMREX_USE_CUB)
template <int BLOCKDIMX>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalAnd_full (int * dest, int source) noexcept
{
    typedef cub::BlockReduce<int,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::TempStorage temp_storage;
    // __syncthreads() prior to writing to shared memory is necessary
    // if this reduction call is occurring multiple times in a kernel,
    // and since we don't know how many times the user is calling it,
    // we do it always to be safe.
    __syncthreads();
    // Compute the block-wide reduction for thread0
    int aggregate = BlockReduce(temp_storage).Reduce(source, amrex::LogicalAnd<int>());
    if (threadIdx.x == 0) { Gpu::Atomic::LogicalAnd(dest, aggregate); }
}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalAnd_full (int * dest, int source) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalAnd<int> >(), 1);
    if (threadIdx.x == 0) { Gpu::Atomic::LogicalAnd(dest, source); }
}

#if defined(AMREX_USE_CUB)
template <int BLOCKDIMX>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalOr_full (int * dest, int source) noexcept
{
    typedef cub::BlockReduce<int,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::TempStorage temp_storage;
    // __syncthreads() prior to writing to shared memory is necessary
    // if this reduction call is occurring multiple times in a kernel,
    // and since we don't know how many times the user is calling it,
    // we do it always to be safe.
    __syncthreads();
    // Compute the block-wide reduction for thread0
    int aggregate = BlockReduce(temp_storage).Reduce(source, amrex::LogicalOr<int>());
    if (threadIdx.x == 0) { Gpu::Atomic::LogicalOr(dest, aggregate); }
}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalOr_full (int * dest, int source) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalOr<int> >(), 0);
    if (threadIdx.x == 0) { Gpu::Atomic::LogicalOr(dest, source); }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceSum (T * dest, T source, Gpu::Handler const& handler) noexcept
{
    if (handler.isFullBlock()) {
#if defined(AMREX_USE_CUB)
        if (blockDim.x == 128) {
            deviceReduceSum_full<128,T>(dest, source);
        } else if (blockDim.x == 256) {
            deviceReduceSum_full<256,T>(dest, source);
        } else if (blockDim.x == 512) {
            deviceReduceSum_full<512,T>(dest, source);
        } else
#endif
        {
            deviceReduceSum_full<T>(dest, source);
        }
    } else {
        Gpu::blockReduce_partial<Gpu::Device::warp_size>
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Plus<T> >(),
             Gpu::AtomicAdd<T>(), handler);
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMin (T * dest, T source, Gpu::Handler const& handler) noexcept
{
    if (handler.isFullBlock()) {
#if defined(AMREX_USE_CUB)
        if (blockDim.x == 128) {
            deviceReduceMin_full<128,T>(dest, source);
        } else if (blockDim.x == 256) {
            deviceReduceMin_full<256,T>(dest, source);
        } else if (blockDim.x == 512) {
            deviceReduceMin_full<512,T>(dest, source);
        } else
#endif
        {
            deviceReduceMin_full<T>(dest, source);
        }
    } else {
        Gpu::blockReduce_partial<Gpu::Device::warp_size>
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Less<T> >(),
             Gpu::AtomicMin<T>(), handler);
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMax (T * dest, T source, Gpu::Handler const& handler) noexcept
{
    if (handler.isFullBlock()) {
#if defined(AMREX_USE_CUB)
        if (blockDim.x == 128) {
            deviceReduceMax_full<128,T>(dest, source);
        } else if (blockDim.x == 256) {
            deviceReduceMax_full<256,T>(dest, source);
        } else if (blockDim.x == 512) {
            deviceReduceMax_full<512,T>(dest, source);
        } else
#endif
        {
            deviceReduceMax_full<T>(dest, source);
        }
    } else {
        Gpu::blockReduce_partial<Gpu::Device::warp_size>
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Greater<T> >(),
             Gpu::AtomicMax<T>(), handler);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalAnd (int * dest, int source, Gpu::Handler const& handler) noexcept
{
    if (handler.isFullBlock()) {
#if defined(AMREX_USE_CUB)
        if (blockDim.x == 128) {
            deviceReduceLogicalAnd_full<128>(dest, source);
        } else if (blockDim.x == 256) {
            deviceReduceLogicalAnd_full<256>(dest, source);
        } else if (blockDim.x == 512) {
            deviceReduceLogicalAnd_full<512>(dest, source);
        } else
#endif
        {
            deviceReduceLogicalAnd_full(dest, source);
        }
    } else {
        Gpu::blockReduce_partial<Gpu::Device::warp_size>
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalAnd<int> >(),
             Gpu::AtomicLogicalAnd<int>(), handler);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalOr (int * dest, int source, Gpu::Handler const& handler) noexcept
{
    if (handler.isFullBlock()) {
#if defined(AMREX_USE_CUB)
        if (blockDim.x == 128) {
            deviceReduceLogicalOr_full<128>(dest, source);
        } else if (blockDim.x == 256) {
            deviceReduceLogicalOr_full<256>(dest, source);
        } else if (blockDim.x == 512) {
            deviceReduceLogicalOr_full<512>(dest, source);
        } else
#endif
        {
            deviceReduceLogicalOr_full(dest, source);
        }
    } else {
        Gpu::blockReduce_partial<Gpu::Device::warp_size>
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalOr<int> >(),
             Gpu::AtomicLogicalOr<int>(), handler);
    }
}

#else

template <typename T>
AMREX_FORCE_INLINE
void deviceReduceSum_full (T * dest, T source) noexcept
{
#ifdef AMREX_USE_OMP
#pragma omp atomic
#endif
    *dest += source;
}

template <typename T>
AMREX_FORCE_INLINE
void deviceReduceSum (T * dest, T source, Gpu::Handler const&) noexcept
{
    deviceReduceSum_full(dest, source);
}

template <typename T>
AMREX_FORCE_INLINE
void deviceReduceMin_full (T * dest, T source) noexcept
{
#ifdef AMREX_USE_OMP
#pragma omp critical (gpureduce_reducemin)
#endif
    *dest = std::min(*dest, source);
}

template <typename T>
AMREX_FORCE_INLINE
void deviceReduceMin (T * dest, T source, Gpu::Handler const&) noexcept
{
    deviceReduceMin_full(dest, source);
}

template <typename T>
AMREX_FORCE_INLINE
void deviceReduceMax_full (T * dest, T source) noexcept
{
#ifdef AMREX_USE_OMP
#pragma omp critical (gpureduce_reducemax)
#endif
    *dest = std::max(*dest, source);
}

template <typename T>
AMREX_FORCE_INLINE
void deviceReduceMax (T * dest, T source, Gpu::Handler const&) noexcept
{
    deviceReduceMax_full(dest, source);
}

AMREX_FORCE_INLINE
void deviceReduceLogicalAnd_full (int * dest, int source) noexcept
{
#ifdef AMREX_USE_OMP
#pragma omp critical (gpureduce_reduceand)
#endif
    *dest = (*dest) && source;
}

AMREX_FORCE_INLINE
void deviceReduceLogicalAnd (int * dest, int source, Gpu::Handler const&) noexcept
{
    deviceReduceLogicalAnd_full(dest, source);
}

AMREX_FORCE_INLINE
void deviceReduceLogicalOr_full (int * dest, int source) noexcept
{
#ifdef AMREX_USE_OMP
#pragma omp critical (gpureduce_reduceor)
#endif
    *dest = (*dest) || source;
}

AMREX_FORCE_INLINE
void deviceReduceLogicalOr (int * dest, int source, Gpu::Handler const&) noexcept
{
    deviceReduceLogicalOr_full(dest, source);
}

#endif
}}

#endif
