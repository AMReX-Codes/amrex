#ifndef AMREX_CUDA_REDUCE_H_
#define AMREX_CUDA_REDUCE_H_

#include <AMReX_GpuQualifiers.H>

namespace amrex {
namespace Cuda {
#ifdef AMREX_USE_CUDA

// Based on https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf by Mark Harris

template <unsigned int blockSize, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void warpxReduceSum (volatile T* data, int tid)
{
    if (blockSize >= 64) {
        data[tid] += data[tid + 32];
#if __CUDA_ARCH__ >= 700
        __syncwarp();
#endif
    }
    if (blockSize >= 32) {
        data[tid] += data[tid + 16];
#if __CUDA_ARCH__ >= 700
        __syncwarp();
#endif
    }
    if (blockSize >= 16) {
        data[tid] += data[tid + 8];
#if __CUDA_ARCH__ >= 700
        __syncwarp();
#endif
    }
    if (blockSize >=  8) {
        data[tid] += data[tid + 4];
#if __CUDA_ARCH__ >= 700
        __syncwarp();
#endif
    }
    if (blockSize >=  4) {
        data[tid] += data[tid + 2];
#if __CUDA_ARCH__ >= 700
        __syncwarp();
#endif
    }
    if (blockSize >=  2) {
        data[tid] += data[tid + 1];
#if __CUDA_ARCH__ >= 700
        __syncwarp();
#endif
    }
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
    if (tid < 32) warpxReduceSum<blockSize>(data, tid);
    if (tid == 0) sum = data[0];
}

#endif
}
}

#endif
