#ifndef AMREX_CUDA_REDUCE_H_
#define AMREX_CUDA_REDUCE_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Utility.H>

namespace amrex {
namespace Cuda {
#ifdef AMREX_USE_CUDA

// Based on https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf by Mark Harris

// sum

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceSum_lt7 (volatile T* data, int tid) noexcept
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
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceSum_ge7 (T* data, int tid) noexcept
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
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceSum (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    warpReduceSum_ge7<blockSize>(data, tid);
#else
    warpReduceSum_lt7<blockSize>(data, tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
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
    if (blockSize >= 128) { if (tid <  64) { data[tid] += data[tid+ 64]; } __syncthreads(); }
    if (tid < 32) warpReduceSum<blockSize>(data, tid);
    if (tid == 0) sum = data[0];
}

// min

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceMin_lt7 (volatile T* data, int tid) noexcept
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
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceMin_ge7 (T* data, int tid) noexcept
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
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceMin (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    warpReduceMin_ge7<blockSize>(data, tid);
#else
    warpReduceMin_lt7<blockSize>(data, tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
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
    if (blockSize >= 128) { if (tid <  64) { data[tid] = amrex::min(data[tid],data[tid+ 64]); } __syncthreads(); }
    if (tid < 32) warpReduceMin<blockSize>(data, tid);
    if (tid == 0) dmin = data[0];
}

// max

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceMax_lt7 (volatile T* data, int tid) noexcept
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
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceMax_ge7 (T* data, int tid) noexcept
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
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceMax (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    warpReduceMax_ge7<blockSize>(data, tid);
#else
    warpReduceMax_lt7<blockSize>(data, tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
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
    if (blockSize >= 128) { if (tid <  64) { data[tid] = amrex::max(data[tid],data[tid+ 64]); } __syncthreads(); }
    if (tid < 32) warpReduceMax<blockSize>(data, tid);
    if (tid == 0) dmax = data[0];
}

// and

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceAnd_lt7 (volatile T* data, int tid) noexcept
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
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceAnd_ge7 (T* data, int tid) noexcept
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
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceAnd (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    warpReduceAnd_ge7<blockSize>(data, tid);
#else
    warpReduceAnd_lt7<blockSize>(data, tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
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
    if (blockSize >= 128) { if (tid <  64) { data[tid] = data[tid] && data[tid+ 64]; } __syncthreads(); }
    if (tid < 32) warpReduceAnd<blockSize>(data, tid);
    if (tid == 0) r = data[0];
}

// or

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceOr_lt7 (volatile T* data, int tid) noexcept
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
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceOr_ge7 (T* data, int tid) noexcept
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
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceOr (T* data, int tid) noexcept
{
#if __CUDA_ARCH__ >= 700
    warpReduceOr_ge7<blockSize>(data, tid);
#else
    warpReduceOr_lt7<blockSize>(data, tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
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
    if (blockSize >= 128) { if (tid <  64) { data[tid] = data[tid] || data[tid+ 64]; } __syncthreads(); }
    if (tid < 32) warpReduceOr<blockSize>(data, tid);
    if (tid == 0) r = data[0];
}

#endif
}
}

#endif
