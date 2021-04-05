#ifndef AMREX_GPU_REDUCE_H_
#define AMREX_GPU_REDUCE_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuTypes.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuUtility.H>
#include <AMReX_Functional.H>

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
    if (lane == 0) shared[wid] = x;
    h.item->barrier(sycl::access::fence_space::local_space);
    bool b = (tid == 0) || (tid < numwarps);
    x =  b ? shared[lane] : x0;
    if (wid == 0) x = warp_reduce(x, sg);
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
       if (sg.get_local_id()[0] == 0) atomic_op(dest, x);
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
    if (h.item->get_local_id(0) == 0) Gpu::Atomic::AddNoRet(dest, source);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMin_full (T * dest, T source, Gpu::Handler const& h) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Less<T> >(), source, h);
    if (h.item->get_local_id(0) == 0) Gpu::Atomic::Min(dest, source);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMax_full (T * dest, T source, Gpu::Handler const& h) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Greater<T> >(), source, h);
    if (h.item->get_local_id(0) == 0) Gpu::Atomic::Max(dest, source);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalAnd_full (int * dest, int source, Gpu::Handler const& h) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalAnd<int> >(), 1, h);
    if (h.item->get_local_id(0) == 0) Gpu::Atomic::LogicalAnd(dest, source);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalOr_full (int * dest, int source, Gpu::Handler const& h) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalOr<int> >(), 0, h);
    if (h.item->get_local_id(0) == 0) Gpu::Atomic::LogicalOr(dest, source);
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
    if (lane == 0) shared[wid] = x;
    __syncthreads();
    bool b = (threadIdx.x == 0) || (threadIdx.x < blockDim.x / warpSize);
    x =  b ? shared[lane] : x0;
    if (wid == 0) x = warp_reduce(x);
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
        if (threadIdx.x % warpSize == 0) atomic_op(dest, x);
    } else {
        atomic_op(dest,x);
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceSum_full (T * dest, T source) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Plus<T> >(), (T)0);
    if (threadIdx.x == 0) Gpu::Atomic::AddNoRet(dest, source);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMin_full (T * dest, T source) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Less<T> >(), source);
    if (threadIdx.x == 0) Gpu::Atomic::Min(dest, source);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMax_full (T * dest, T source) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,T,amrex::Greater<T> >(), source);
    if (threadIdx.x == 0) Gpu::Atomic::Max(dest, source);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalAnd_full (int * dest, int source) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalAnd<int> >(), 1);
    if (threadIdx.x == 0) Gpu::Atomic::LogicalAnd(dest, source);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceLogicalOr_full (int * dest, int source) noexcept
{
    source = Gpu::blockReduce<Gpu::Device::warp_size>
        (source, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::LogicalOr<int> >(), 0);
    if (threadIdx.x == 0) Gpu::Atomic::LogicalOr(dest, source);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceSum (T * dest, T source, Gpu::Handler const& handler) noexcept
{
    if (handler.isFullBlock()) {
        deviceReduceSum_full(dest, source);
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
        deviceReduceMin_full(dest, source);
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
        deviceReduceMax_full(dest, source);
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
        deviceReduceLogicalAnd_full(dest, source);
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
        deviceReduceLogicalOr_full(dest, source);
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

namespace amrex {
namespace Gpu {
#if defined(AMREX_USE_GPU) && !defined(AMREX_USE_DPCPP)

// Based on https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf by Mark Harris

// sum

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void amdWarpReduceSum (volatile T* data, int tid) noexcept
{
#ifdef __HIP_DEVICE_COMPILE__
    if (blockSize >= 128) data[tid] += data[tid + 64];
    if (blockSize >= 64) data[tid] += data[tid + 32];
    if (blockSize >= 32) data[tid] += data[tid + 16];
    if (blockSize >= 16) data[tid] += data[tid + 8];
    if (blockSize >=  8) data[tid] += data[tid + 4];
    if (blockSize >=  4) data[tid] += data[tid + 2];
    if (blockSize >=  2) data[tid] += data[tid + 1];
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceSum_lt7 (volatile T* data, int tid) noexcept
{
#if __CUDA_ARCH__ < 700
    if (blockSize >= 64) data[tid] += data[tid + 32];
    if (blockSize >= 32) data[tid] += data[tid + 16];
    if (blockSize >= 16) data[tid] += data[tid + 8];
    if (blockSize >=  8) data[tid] += data[tid + 4];
    if (blockSize >=  4) data[tid] += data[tid + 2];
    if (blockSize >=  2) data[tid] += data[tid + 1];
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceSum_ge7 (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    if (blockSize >= 64) { if (tid < 32) { data[tid] += data[tid + 32]; } __syncwarp(); }
    if (blockSize >= 32) { if (tid < 16) { data[tid] += data[tid + 16]; } __syncwarp(); }
    if (blockSize >= 16) { if (tid <  8) { data[tid] += data[tid +  8]; } __syncwarp(); }
    if (blockSize >=  8) { if (tid <  4) { data[tid] += data[tid +  4]; } __syncwarp(); }
    if (blockSize >=  4) { if (tid <  2) { data[tid] += data[tid +  2]; } __syncwarp(); }
    if (blockSize >=  2) { if (tid <  1) { data[tid] += data[tid +  1]; } __syncwarp(); }
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceSum (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    cudaWarpReduceSum_ge7<blockSize>(data, tid);
#else
    cudaWarpReduceSum_lt7<blockSize>(data, tid);
#endif
}

template <unsigned int blockSize, int warpSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void blockReduceSum (T* data, T& sum) noexcept
{
    int tid = threadIdx.x;
    if (blockSize >= 1024) {
        if (tid < 512) {
            for (int n = tid+512; n < blockSize; n += 512) {
                data[tid] += data[n];
            }
        }
        __syncthreads();
    }
    if (blockSize >= 512) { if (tid < 256) { data[tid] += data[tid+256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { data[tid] += data[tid+128]; } __syncthreads(); }
    if (warpSize >= 64) {
        if (tid < 64) amdWarpReduceSum<blockSize>(data, tid);
    } else {
        if (blockSize >= 128) { if (tid <  64) { data[tid] += data[tid+ 64]; } __syncthreads(); }
        if (tid < 32) cudaWarpReduceSum<blockSize>(data, tid);
    }
    if (tid == 0) sum = data[0];
}

// min

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void amdWarpReduceMin (volatile T* data, int tid) noexcept
{
#ifdef __HIP_DEVICE_COMPILE__
    if (blockSize >= 128) data[tid] = amrex::min(data[tid],data[tid + 64]);
    if (blockSize >= 64) data[tid] = amrex::min(data[tid],data[tid + 32]);
    if (blockSize >= 32) data[tid] = amrex::min(data[tid],data[tid + 16]);
    if (blockSize >= 16) data[tid] = amrex::min(data[tid],data[tid +  8]);
    if (blockSize >=  8) data[tid] = amrex::min(data[tid],data[tid +  4]);
    if (blockSize >=  4) data[tid] = amrex::min(data[tid],data[tid +  2]);
    if (blockSize >=  2) data[tid] = amrex::min(data[tid],data[tid +  1]);
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceMin_lt7 (volatile T* data, int tid) noexcept
{
#if __CUDA_ARCH__ < 700
    if (blockSize >= 64) data[tid] = amrex::min(data[tid],data[tid + 32]);
    if (blockSize >= 32) data[tid] = amrex::min(data[tid],data[tid + 16]);
    if (blockSize >= 16) data[tid] = amrex::min(data[tid],data[tid +  8]);
    if (blockSize >=  8) data[tid] = amrex::min(data[tid],data[tid +  4]);
    if (blockSize >=  4) data[tid] = amrex::min(data[tid],data[tid +  2]);
    if (blockSize >=  2) data[tid] = amrex::min(data[tid],data[tid +  1]);
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceMin_ge7 (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    if (blockSize >= 64) { if (tid < 32) { data[tid] = amrex::min(data[tid],data[tid + 32]); } __syncwarp(); }
    if (blockSize >= 32) { if (tid < 16) { data[tid] = amrex::min(data[tid],data[tid + 16]); } __syncwarp(); }
    if (blockSize >= 16) { if (tid <  8) { data[tid] = amrex::min(data[tid],data[tid +  8]); } __syncwarp(); }
    if (blockSize >=  8) { if (tid <  4) { data[tid] = amrex::min(data[tid],data[tid +  4]); } __syncwarp(); }
    if (blockSize >=  4) { if (tid <  2) { data[tid] = amrex::min(data[tid],data[tid +  2]); } __syncwarp(); }
    if (blockSize >=  2) { if (tid <  1) { data[tid] = amrex::min(data[tid],data[tid +  1]); } __syncwarp(); }
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceMin (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    cudaWarpReduceMin_ge7<blockSize>(data, tid);
#else
    cudaWarpReduceMin_lt7<blockSize>(data, tid);
#endif
}

template <unsigned int blockSize, int warpSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void blockReduceMin (T* data, T& dmin) noexcept
{
    int tid = threadIdx.x;
    if (blockSize >= 1024) {
        if (tid < 512) {
            for (int n = tid+512; n < blockSize; n += 512) {
                data[tid] = amrex::min(data[tid],data[n]);
            }
        }
        __syncthreads();
    }
    if (blockSize >= 512) { if (tid < 256) { data[tid] = amrex::min(data[tid],data[tid+256]); } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { data[tid] = amrex::min(data[tid],data[tid+128]); } __syncthreads(); }
    if (warpSize >= 64) {
        if (tid < 64) amdWarpReduceMin<blockSize>(data, tid);
    } else {
        if (blockSize >= 128) { if (tid <  64) { data[tid] = amrex::min(data[tid],data[tid+ 64]); } __syncthreads(); }
        if (tid < 32) cudaWarpReduceMin<blockSize>(data, tid);
    }
    if (tid == 0) dmin = data[0];
}

// max

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void amdWarpReduceMax (volatile T* data, int tid) noexcept
{
#ifdef __HIP_DEVICE_COMPILE__
    if (blockSize >= 128) data[tid] = amrex::max(data[tid],data[tid + 64]);
    if (blockSize >= 64) data[tid] = amrex::max(data[tid],data[tid + 32]);
    if (blockSize >= 32) data[tid] = amrex::max(data[tid],data[tid + 16]);
    if (blockSize >= 16) data[tid] = amrex::max(data[tid],data[tid +  8]);
    if (blockSize >=  8) data[tid] = amrex::max(data[tid],data[tid +  4]);
    if (blockSize >=  4) data[tid] = amrex::max(data[tid],data[tid +  2]);
    if (blockSize >=  2) data[tid] = amrex::max(data[tid],data[tid +  1]);
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceMax_lt7 (volatile T* data, int tid) noexcept
{
#if __CUDA_ARCH__ < 700
    if (blockSize >= 64) data[tid] = amrex::max(data[tid],data[tid + 32]);
    if (blockSize >= 32) data[tid] = amrex::max(data[tid],data[tid + 16]);
    if (blockSize >= 16) data[tid] = amrex::max(data[tid],data[tid +  8]);
    if (blockSize >=  8) data[tid] = amrex::max(data[tid],data[tid +  4]);
    if (blockSize >=  4) data[tid] = amrex::max(data[tid],data[tid +  2]);
    if (blockSize >=  2) data[tid] = amrex::max(data[tid],data[tid +  1]);
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceMax_ge7 (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    if (blockSize >= 64) { if (tid < 32) { data[tid] = amrex::max(data[tid],data[tid + 32]); } __syncwarp(); }
    if (blockSize >= 32) { if (tid < 16) { data[tid] = amrex::max(data[tid],data[tid + 16]); } __syncwarp(); }
    if (blockSize >= 16) { if (tid <  8) { data[tid] = amrex::max(data[tid],data[tid +  8]); } __syncwarp(); }
    if (blockSize >=  8) { if (tid <  4) { data[tid] = amrex::max(data[tid],data[tid +  4]); } __syncwarp(); }
    if (blockSize >=  4) { if (tid <  2) { data[tid] = amrex::max(data[tid],data[tid +  2]); } __syncwarp(); }
    if (blockSize >=  2) { if (tid <  1) { data[tid] = amrex::max(data[tid],data[tid +  1]); } __syncwarp(); }
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceMax (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    cudaWarpReduceMax_ge7<blockSize>(data, tid);
#else
    cudaWarpReduceMax_lt7<blockSize>(data, tid);
#endif
}

template <unsigned int blockSize, int warpSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void blockReduceMax (T* data, T& dmax) noexcept
{
    int tid = threadIdx.x;
    if (blockSize >= 1024) {
        if (tid < 512) {
            for (int n = tid+512; n < blockSize; n += 512) {
                data[tid] = amrex::max(data[tid],data[n]);
            }
        }
        __syncthreads();
    }
    if (blockSize >= 512) { if (tid < 256) { data[tid] = amrex::max(data[tid],data[tid+256]); } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { data[tid] = amrex::max(data[tid],data[tid+128]); } __syncthreads(); }
    if (warpSize >= 64) {
        if (tid < 64) amdWarpReduceMax<blockSize>(data, tid);
    } else {
        if (blockSize >= 128) { if (tid <  64) { data[tid] = amrex::max(data[tid],data[tid+ 64]); } __syncthreads(); }
        if (tid < 32) cudaWarpReduceMax<blockSize>(data, tid);
    }
    if (tid == 0) dmax = data[0];
}

// and

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void amdWarpReduceAnd (volatile T* data, int tid) noexcept
{
#ifdef __HIP_DEVICE_COMPILE__
    if (blockSize >= 128) data[tid] = data[tid] && data[tid + 64];
    if (blockSize >= 64) data[tid] = data[tid] && data[tid + 32];
    if (blockSize >= 32) data[tid] = data[tid] && data[tid + 16];
    if (blockSize >= 16) data[tid] = data[tid] && data[tid +  8];
    if (blockSize >=  8) data[tid] = data[tid] && data[tid +  4];
    if (blockSize >=  4) data[tid] = data[tid] && data[tid +  2];
    if (blockSize >=  2) data[tid] = data[tid] && data[tid +  1];
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceAnd_lt7 (volatile T* data, int tid) noexcept
{
#if __CUDA_ARCH__ < 700
    if (blockSize >= 64) data[tid] = data[tid] && data[tid + 32];
    if (blockSize >= 32) data[tid] = data[tid] && data[tid + 16];
    if (blockSize >= 16) data[tid] = data[tid] && data[tid +  8];
    if (blockSize >=  8) data[tid] = data[tid] && data[tid +  4];
    if (blockSize >=  4) data[tid] = data[tid] && data[tid +  2];
    if (blockSize >=  2) data[tid] = data[tid] && data[tid +  1];
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceAnd_ge7 (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    if (blockSize >= 64) { if (tid < 32) { data[tid] = data[tid] && data[tid + 32]; } __syncwarp(); }
    if (blockSize >= 32) { if (tid < 16) { data[tid] = data[tid] && data[tid + 16]; } __syncwarp(); }
    if (blockSize >= 16) { if (tid <  8) { data[tid] = data[tid] && data[tid +  8]; } __syncwarp(); }
    if (blockSize >=  8) { if (tid <  4) { data[tid] = data[tid] && data[tid +  4]; } __syncwarp(); }
    if (blockSize >=  4) { if (tid <  2) { data[tid] = data[tid] && data[tid +  2]; } __syncwarp(); }
    if (blockSize >=  2) { if (tid <  1) { data[tid] = data[tid] && data[tid +  1]; } __syncwarp(); }
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceAnd (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    cudaWarpReduceAnd_ge7<blockSize>(data, tid);
#else
    cudaWarpReduceAnd_lt7<blockSize>(data, tid);
#endif
}

template <unsigned int blockSize, int warpSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void blockReduceAnd (T* data, T& r) noexcept
{
    int tid = threadIdx.x;
    if (blockSize >= 1024) {
        if (tid < 512) {
            for (int n = tid+512; n < blockSize; n += 512) {
                data[tid] = data[tid] && data[n];
            }
        }
        __syncthreads();
    }
    if (blockSize >= 512) { if (tid < 256) { data[tid] = data[tid] && data[tid+256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { data[tid] = data[tid] && data[tid+128]; } __syncthreads(); }
    if (warpSize >= 64) {
        if (tid < 64) amdWarpReduceAnd<blockSize>(data, tid);
    } else {
        if (blockSize >= 128) { if (tid <  64) { data[tid] = data[tid] && data[tid+ 64]; } __syncthreads(); }
        if (tid < 32) cudaWarpReduceAnd<blockSize>(data, tid);
    }
    if (tid == 0) r = data[0];
}

// or

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void amdWarpReduceOr (volatile T* data, int tid) noexcept
{
#ifdef __HIP_DEVICE_COMPILE__
    if (blockSize >= 128) data[tid] = data[tid] || data[tid + 64];
    if (blockSize >= 64) data[tid] = data[tid] || data[tid + 32];
    if (blockSize >= 32) data[tid] = data[tid] || data[tid + 16];
    if (blockSize >= 16) data[tid] = data[tid] || data[tid +  8];
    if (blockSize >=  8) data[tid] = data[tid] || data[tid +  4];
    if (blockSize >=  4) data[tid] = data[tid] || data[tid +  2];
    if (blockSize >=  2) data[tid] = data[tid] || data[tid +  1];
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceOr_lt7 (volatile T* data, int tid) noexcept
{
#if __CUDA_ARCH__ < 700
    if (blockSize >= 64) data[tid] = data[tid] || data[tid + 32];
    if (blockSize >= 32) data[tid] = data[tid] || data[tid + 16];
    if (blockSize >= 16) data[tid] = data[tid] || data[tid +  8];
    if (blockSize >=  8) data[tid] = data[tid] || data[tid +  4];
    if (blockSize >=  4) data[tid] = data[tid] || data[tid +  2];
    if (blockSize >=  2) data[tid] = data[tid] || data[tid +  1];
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceOr_ge7 (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    if (blockSize >= 64) { if (tid < 32) { data[tid] = data[tid] || data[tid + 32]; } __syncwarp(); }
    if (blockSize >= 32) { if (tid < 16) { data[tid] = data[tid] || data[tid + 16]; } __syncwarp(); }
    if (blockSize >= 16) { if (tid <  8) { data[tid] = data[tid] || data[tid +  8]; } __syncwarp(); }
    if (blockSize >=  8) { if (tid <  4) { data[tid] = data[tid] || data[tid +  4]; } __syncwarp(); }
    if (blockSize >=  4) { if (tid <  2) { data[tid] = data[tid] || data[tid +  2]; } __syncwarp(); }
    if (blockSize >=  2) { if (tid <  1) { data[tid] = data[tid] || data[tid +  1]; } __syncwarp(); }
#else
    amrex::ignore_unused(data,tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cudaWarpReduceOr (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    cudaWarpReduceOr_ge7<blockSize>(data, tid);
#else
    cudaWarpReduceOr_lt7<blockSize>(data, tid);
#endif
}

template <unsigned int blockSize, int warpSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void blockReduceOr (T* data, T& r) noexcept
{
    int tid = threadIdx.x;
    if (blockSize >= 1024) {
        if (tid < 512) {
            for (int n = tid+512; n < blockSize; n += 512) {
                data[tid] = data[tid] || data[n];
            }
        }
        __syncthreads();
    }
    if (blockSize >= 512) { if (tid < 256) { data[tid] = data[tid] || data[tid+256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { data[tid] = data[tid] || data[tid+128]; } __syncthreads(); }
    if (warpSize >= 64) {
        if (tid < 64) amdWarpReduceOr<blockSize>(data, tid);
    } else {
        if (blockSize >= 128) { if (tid <  64) { data[tid] = data[tid] || data[tid+ 64]; } __syncthreads(); }
        if (tid < 32) cudaWarpReduceOr<blockSize>(data, tid);
    }
    if (tid == 0) r = data[0];
}

#endif
}
}

#endif
