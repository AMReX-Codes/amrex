#ifndef AMREX_CUDA_REDUCE_H_
#define AMREX_CUDA_REDUCE_H_

#include <AMReX_GpuQualifiers.H>

namespace amrex {
namespace Cuda {
#ifdef AMREX_USE_CUDA

// Based on https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf by Mark Harris

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void warpReduceSum_lt7 (volatile T* data, int tid)
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
void warpReduceSum_ge7 (T* data, int tid)
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
void warpReduceSum (T* data, int tid)
{
#if __CUDA_ARCH__ >= 700
    warpReduceSum_ge7<blockSize>(data, tid);
#else
    warpReduceSum_lt7<blockSize>(data, tid);
#endif
}

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE
void blockReduceSum (T* data, T& sum)
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

#endif
}
}

#endif
