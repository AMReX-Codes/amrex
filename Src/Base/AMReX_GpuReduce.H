#ifndef AMREX_GPU_REDUCE_H_
#define AMREX_GPU_REDUCE_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuTypes.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuUtility.H>
#include <AMReX_Functional.H>
#include <AMReX_TypeTraits.H>

#if !defined(AMREX_USE_CUB) && defined(AMREX_USE_CUDA) && defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 11)
#define AMREX_USE_CUB 1
#endif

#if defined(AMREX_USE_CUB)
#include <cub/cub.cuh>
#elif defined(AMREX_USE_HIP)
#include <rocprim/rocprim.hpp>
#endif

//
// Public interface
//
namespace amrex::Gpu {
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
}

//
// Reduce functions based on _shfl_down_sync
//

namespace amrex::Gpu {

#ifdef AMREX_USE_SYCL

template <int warpSize, typename T, typename F>
struct warpReduce
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T operator() (T x, sycl::sub_group const& sg) const noexcept
    {
        for (int offset = warpSize/2; offset > 0; offset /= 2) {
            T y = sycl::shift_group_left(sg, x, offset);
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
    sycl::sub_group const& sg = h.item->get_sub_group();
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
   sycl::sub_group const& sg = handler.item->get_sub_group();
   int wid = sg.get_group_id()[0];
   if ((wid+1)*warpSize <= handler.numActiveThreads) {
       x = warp_reduce(x, sg); // full warp
       if (sg.get_local_id()[0] == 0) { atomic_op(dest, x); }
   } else {
       atomic_op(dest, x);
   }
}

// block-wide reduction for thread0
template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduceSum (T source, Gpu::Handler const& h) noexcept
{
    return Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Plus<T> >(), (T)0, h);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceSum_full (T * dest, T source, Gpu::Handler const& h) noexcept
{
    T aggregate = blockReduceSum(source, h);
    if (h.item->get_local_id(0) == 0) { Gpu::Atomic::AddNoRet(dest, aggregate); }
}

// block-wide reduction for thread0
template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduceMin (T source, Gpu::Handler const& h) noexcept
{
    return Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Minimum<T> >(), source, h);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMin_full (T * dest, T source, Gpu::Handler const& h) noexcept
{
    T aggregate = blockReduceMin(source, h);
    if (h.item->get_local_id(0) == 0) { Gpu::Atomic::Min(dest, aggregate); }
}

// block-wide reduction for thread0
template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduceMax (T source, Gpu::Handler const& h) noexcept
{
    return Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Maximum<T> >(), source, h);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMax_full (T * dest, T source, Gpu::Handler const& h) noexcept
{
    T aggregate = blockReduceMax(source, h);
    if (h.item->get_local_id(0) == 0) { Gpu::Atomic::Max(dest, aggregate); }
}

// block-wide reduction for thread0
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int blockReduceLogicalAnd (int source, Gpu::Handler const& h) noexcept
{
    return Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalAnd<int> >(), 1, h);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalAnd_full (int * dest, int source, Gpu::Handler const& h) noexcept
{
    int aggregate = blockReduceLogicalAnd(source, h);
    if (h.item->get_local_id(0) == 0) { Gpu::Atomic::LogicalAnd(dest, aggregate); }
}

// block-wide reduction for thread0
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int blockReduceLogicalOr (int source, Gpu::Handler const& h) noexcept
{
    return Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalOr<int> >(), 0, h);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalOr_full (int * dest, int source, Gpu::Handler const& h) noexcept
{
    int aggregate = blockReduceLogicalOr(source, h);
    if (h.item->get_local_id(0) == 0) { Gpu::Atomic::LogicalOr(dest, aggregate); }
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
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Minimum<T> >(),
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
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Maximum<T> >(),
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

namespace detail {

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T shuffle_down (T x, int offset) noexcept
{
    return AMREX_HIP_OR_CUDA(__shfl_down(x, offset),
                             __shfl_down_sync(0xffffffff, x, offset));
}

// If other sizeof is needed, we can implement it later.
template <class T, std::enable_if_t<sizeof(T)%sizeof(unsigned int) == 0, int> = 0>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T multi_shuffle_down (T x, int offset) noexcept
{
    constexpr int nwords = (sizeof(T) + sizeof(unsigned int) - 1) / sizeof(unsigned int);
    T y;
    auto py = reinterpret_cast<unsigned int*>(&y);
    auto px = reinterpret_cast<unsigned int*>(&x);
    for (int i = 0; i < nwords; ++i) {
        py[i] = shuffle_down(px[i],offset);
    }
    return y;
}

}

template <int warpSize, typename T, typename F>
struct warpReduce
{
    // Not all arithmetic types can be taken by shuffle_down, but it's good enough.
    template <class U=T, std::enable_if_t<std::is_arithmetic<U>::value,int> = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T operator() (T x) const noexcept
    {
        for (int offset = warpSize/2; offset > 0; offset /= 2) {
            T y = detail::shuffle_down(x, offset);
            x = F()(x,y);
        }
        return x;
    }

    template <class U=T, std::enable_if_t<!std::is_arithmetic<U>::value,int> = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T operator() (T x) const noexcept
    {
        for (int offset = warpSize/2; offset > 0; offset /= 2) {
            T y = detail::multi_shuffle_down(x, offset);
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

// Compute the block-wide reduction for thread0
template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduceSum (T source) noexcept
{
    return Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Plus<T> >(), (T)0);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceSum_full (T * dest, T source) noexcept
{
    T aggregate = blockReduceSum(source);
    if (threadIdx.x == 0) { Gpu::Atomic::AddNoRet(dest, aggregate); }
}

// Compute the block-wide reduction for thread0
template <int BLOCKDIMX, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduceSum (T source) noexcept
{
#if defined(AMREX_USE_CUB)
    typedef cub::BlockReduce<T,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::TempStorage temp_storage;
    // __syncthreads() prior to writing to shared memory is necessary
    // if this reduction call is occurring multiple times in a kernel,
    // and since we don't know how many times the user is calling it,
    // we do it always to be safe.
    __syncthreads();
    return BlockReduce(temp_storage).Sum(source);
#elif defined(AMREX_USE_HIP)
    typedef rocprim::block_reduce<T,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::storage_type temp_storage;
    __syncthreads();
    BlockReduce().reduce(source, source, temp_storage, rocprim::plus<T>());
    return source;
#else
    return blockReduceSum(source);
#endif
}

template <int BLOCKDIMX, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceSum_full (T * dest, T source) noexcept
{
    T aggregate = blockReduceSum<BLOCKDIMX>(source);
    if (threadIdx.x == 0) { Gpu::Atomic::AddNoRet(dest, aggregate); }
}

// Compute the block-wide reduction for thread0
template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduceMin (T source) noexcept
{
    return Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Minimum<T> >(), source);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMin_full (T * dest, T source) noexcept
{
    T aggregate = blockReduceMin(source);
    if (threadIdx.x == 0) { Gpu::Atomic::Min(dest, aggregate); }
}

// Compute the block-wide reduction for thread0
template <int BLOCKDIMX, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduceMin (T source) noexcept
{
#if defined(AMREX_USE_CUB)
    typedef cub::BlockReduce<T,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::TempStorage temp_storage;
    // __syncthreads() prior to writing to shared memory is necessary
    // if this reduction call is occurring multiple times in a kernel,
    // and since we don't know how many times the user is calling it,
    // we do it always to be safe.
    __syncthreads();
    return BlockReduce(temp_storage).Reduce(source, cub::Min());
#elif defined(AMREX_USE_HIP)
    typedef rocprim::block_reduce<T,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::storage_type temp_storage;
    __syncthreads();
    BlockReduce().reduce(source, source, temp_storage, rocprim::minimum<T>());
    return source;
#else
    return blockReduceMin(source);
#endif
}

template <int BLOCKDIMX, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMin_full (T * dest, T source) noexcept
{
    T aggregate = blockReduceMin<BLOCKDIMX>(source);
    if (threadIdx.x == 0) { Gpu::Atomic::Min(dest, aggregate); }
}

// Compute the block-wide reduction for thread0
template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduceMax (T source) noexcept
{
    return Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Maximum<T> >(), source);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMax_full (T * dest, T source) noexcept
{
    T aggregate = blockReduceMax(source);
    if (threadIdx.x == 0) { Gpu::Atomic::Max(dest, aggregate); }
}

// Compute the block-wide reduction for thread0
template <int BLOCKDIMX, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduceMax (T source) noexcept
{
#if defined(AMREX_USE_CUB)
    typedef cub::BlockReduce<T,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::TempStorage temp_storage;
    // __syncthreads() prior to writing to shared memory is necessary
    // if this reduction call is occurring multiple times in a kernel,
    // and since we don't know how many times the user is calling it,
    // we do it always to be safe.
    __syncthreads();
    return BlockReduce(temp_storage).Reduce(source, cub::Max());
#elif defined(AMREX_USE_HIP)
    typedef rocprim::block_reduce<T,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::storage_type temp_storage;
    __syncthreads();
    BlockReduce().reduce(source, source, temp_storage, rocprim::maximum<T>());
    return source;
#else
    return blockReduceMax(source);
#endif
}

template <int BLOCKDIMX, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMax_full (T * dest, T source) noexcept
{
    T aggregate = blockReduceMax<BLOCKDIMX>(source);
    if (threadIdx.x == 0) { Gpu::Atomic::Max(dest, aggregate); }
}

// \cond CODEGEN

// Compute the block-wide reduction for thread0
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int blockReduceLogicalAnd (int source) noexcept
{
    return Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalAnd<int> >(), 1);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalAnd_full (int * dest, int source) noexcept
{
    int aggregate = blockReduceLogicalAnd(source);
    if (threadIdx.x == 0) { Gpu::Atomic::LogicalAnd(dest, aggregate); }
}

// Compute the block-wide reduction for thread0
template <int BLOCKDIMX>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int blockReduceLogicalAnd (int source) noexcept
{
#if defined(AMREX_USE_CUB)
    typedef cub::BlockReduce<int,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::TempStorage temp_storage;
    // __syncthreads() prior to writing to shared memory is necessary
    // if this reduction call is occurring multiple times in a kernel,
    // and since we don't know how many times the user is calling it,
    // we do it always to be safe.
    __syncthreads();
    return BlockReduce(temp_storage).Reduce(source, amrex::LogicalAnd<int>());
#elif defined(AMREX_USE_HIP)
    typedef rocprim::block_reduce<int,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::storage_type temp_storage;
    __syncthreads();
    BlockReduce().reduce(source, source, temp_storage, amrex::LogicalAnd<int>());
    return source;
#else
    return blockReduceLogicalAnd(source);
#endif
}

template <int BLOCKDIMX>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalAnd_full (int * dest, int source) noexcept
{
    int aggregate = blockReduceLogicalAnd<BLOCKDIMX>(source);
    if (threadIdx.x == 0) { Gpu::Atomic::LogicalAnd(dest, aggregate); }
}

// Compute the block-wide reduction for thread0
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int blockReduceLogicalOr (int source) noexcept
{
    return Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalOr<int> >(), 0);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalOr_full (int * dest, int source) noexcept
{
    int aggregate = blockReduceLogicalOr(source);
    if (threadIdx.x == 0) { Gpu::Atomic::LogicalOr(dest, aggregate); }
}

// Compute the block-wide reduction for thread0
template <int BLOCKDIMX>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int blockReduceLogicalOr (int source) noexcept
{
#if defined(AMREX_USE_CUB)
    typedef cub::BlockReduce<int,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::TempStorage temp_storage;
    // __syncthreads() prior to writing to shared memory is necessary
    // if this reduction call is occurring multiple times in a kernel,
    // and since we don't know how many times the user is calling it,
    // we do it always to be safe.
    __syncthreads();
    return BlockReduce(temp_storage).Reduce(source, amrex::LogicalOr<int>());
#elif defined(AMREX_USE_HIP)
    typedef rocprim::block_reduce<int,BLOCKDIMX> BlockReduce;
    __shared__ typename BlockReduce::storage_type temp_storage;
    __syncthreads();
    BlockReduce().reduce(source, source, temp_storage, amrex::LogicalOr<int>());
    return source;
#else
    return blockReduceLogicalOr(source);
#endif
}

template <int BLOCKDIMX>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalOr_full (int * dest, int source) noexcept
{
    int aggregate = blockReduceLogicalOr<BLOCKDIMX>(source);
    if (threadIdx.x == 0) { Gpu::Atomic::LogicalOr(dest, aggregate); }
}

// \endcond

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceSum (T * dest, T source, Gpu::Handler const& handler) noexcept
{
    if (handler.isFullBlock()) {
        if (blockDim.x == 128) {
            deviceReduceSum_full<128,T>(dest, source);
        } else if (blockDim.x == 256) {
            deviceReduceSum_full<256,T>(dest, source);
        } else {
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
        if (blockDim.x == 128) {
            deviceReduceMin_full<128,T>(dest, source);
        } else if (blockDim.x == 256) {
            deviceReduceMin_full<256,T>(dest, source);
        } else {
            deviceReduceMin_full<T>(dest, source);
        }
    } else {
        Gpu::blockReduce_partial<Gpu::Device::warp_size>
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Minimum<T> >(),
             Gpu::AtomicMin<T>(), handler);
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMax (T * dest, T source, Gpu::Handler const& handler) noexcept
{
    if (handler.isFullBlock()) {
        if (blockDim.x == 128) {
            deviceReduceMax_full<128,T>(dest, source);
        } else if (blockDim.x == 256) {
            deviceReduceMax_full<256,T>(dest, source);
        } else {
            deviceReduceMax_full<T>(dest, source);
        }
    } else {
        Gpu::blockReduce_partial<Gpu::Device::warp_size>
            (dest, source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Maximum<T> >(),
             Gpu::AtomicMax<T>(), handler);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalAnd (int * dest, int source, Gpu::Handler const& handler) noexcept
{
    if (handler.isFullBlock()) {
        if (blockDim.x == 128) {
            deviceReduceLogicalAnd_full<128>(dest, source);
        } else if (blockDim.x == 256) {
            deviceReduceLogicalAnd_full<256>(dest, source);
        } else {
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
        if (blockDim.x == 128) {
            deviceReduceLogicalOr_full<128>(dest, source);
        } else if (blockDim.x == 256) {
            deviceReduceLogicalOr_full<256>(dest, source);
        } else {
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
}

#endif
