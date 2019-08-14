#ifndef AMREX_GPU_REDUCE_H_
#define AMREX_GPU_REDUCE_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuUtility.H>
#include <AMReX_Functional.H>

//
// Reduce functions based on _shfl_down_sync
//

namespace amrex { namespace Gpu {

#ifdef AMREX_USE_GPU

template <int warpSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T warpReduceSum (T x) noexcept
{
    for (int offset = warpSize/2; offset > 0; offset /= 2) {
        x += __shfl_down_sync(0xffffffff, x, offset);
    }
    return x;
}

template <int warpSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduceSum (T x) noexcept
{
    __shared__ T shared[warpSize];
    int lane = threadIdx.x % warpSize;
    int wid = threadIdx.x / warpSize;
    x = warpReduceSum<warpSize,T>(x);
    if (lane == 0) shared[wid] = x;
    __syncthreads();
    bool b = (threadIdx.x < blockDim.x / warpSize);
    x =  b ? shared[lane] : 0;
    if (wid == 0) x = warpReduceSum<warpSize,T>(x);
    return x;
}

template <int warpSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T warpReduceMin (T x) noexcept
{
    for (int offset = warpSize/2; offset > 0; offset /= 2) {
        T y = __shfl_down_sync(0xffffffff, x, offset);
        x = amrex::min(x,y);
    }
    return x;
}

template <int warpSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduceMin (T x) noexcept
{
    __shared__ T shared[warpSize];
    int lane = threadIdx.x % warpSize;
    int wid = threadIdx.x / warpSize;
    x = warpReduceMin<warpSize,T>(x);
    if (lane == 0) shared[wid] = x;
    __syncthreads();
    bool b = (threadIdx.x < blockDim.x / warpSize);
    x =  b ? shared[lane] : x;
    if (wid == 0) x = warpReduceMin<warpSize,T>(x);
    return x;
}

template <int warpSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T warpReduceMax (T x) noexcept
{
    for (int offset = warpSize/2; offset > 0; offset /= 2) {
        T y = __shfl_down_sync(0xffffffff, x, offset);
        x = amrex::max(x,y);
    }
    return x;
}

template <int warpSize, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
T blockReduceMax (T x) noexcept
{
    __shared__ T shared[warpSize];
    int lane = threadIdx.x % warpSize;
    int wid = threadIdx.x / warpSize;
    x = warpReduceMax<warpSize,T>(x);
    if (lane == 0) shared[wid] = x;
    __syncthreads();
    bool b = (threadIdx.x < blockDim.x / warpSize);
    x =  b ? shared[lane] : x;
    if (wid == 0) x = warpReduceMax<warpSize,T>(x);
    return x;
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceSum (T * dest, T source) noexcept
{
    source = Gpu::blockReduceSum<Gpu::Device::warp_size,T>(source);
    if (threadIdx.x == 0) Gpu::Atomic::Add(dest, source);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMin (T * dest, T source) noexcept
{
    source = Gpu::blockReduceMin<Gpu::Device::warp_size,T>(source);
    if (threadIdx.x == 0) Gpu::Atomic::Min(dest, source);
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void deviceReduceMax (T * dest, T source) noexcept
{
    source = Gpu::blockReduceMax<Gpu::Device::warp_size,T>(source);
    if (threadIdx.x == 0) Gpu::Atomic::Max(dest, source);
}

#else

template <typename T>
AMREX_FORCE_INLINE
void deviceReduceSum (T * dest, T source) noexcept
{
#ifdef _OPENMP
#pragma omp atomic
#endif
    *dest += source;
}

template <typename T>
AMREX_FORCE_INLINE
void deviceReduceMin (T * dest, T source) noexcept
{
#ifdef _OPENMP
#pragma omp critical (gpureduce_reducemin)
#endif
    *dest = std::min(*dest, source);
}

template <typename T>
AMREX_FORCE_INLINE
void deviceReduceMax (T * dest, T source) noexcept
{
#ifdef _OPENMP
#pragma omp critical (gpureduce_reducemax)
#endif
    *dest = std::max(*dest, source);
}

#endif
}}

namespace amrex {
namespace Gpu {
#ifdef AMREX_USE_GPU

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
