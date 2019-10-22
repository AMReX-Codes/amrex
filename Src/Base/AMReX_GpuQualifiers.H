#ifndef AMREX_GPU_QUALIFIERS_H_
#define AMREX_GPU_QUALIFIERS_H_

#if defined(AMREX_USE_GPU)

#if defined(AMREX_USE_HIP)
#include <hip/hip_runtime.h>
#endif

#define AMREX_GPU_HOST __host__
#define AMREX_GPU_DEVICE __device__
#define AMREX_GPU_GLOBAL __global__
#define AMREX_GPU_HOST_DEVICE __host__ __device__

#if defined(AMREX_USE_HIP) && !defined(AMREX_USE_CUDA)
#define AMREX_GPU_DEVICE_MANAGED __device__
#else
#define AMREX_GPU_DEVICE_MANAGED __device__ __managed__
#endif

// Put in this header because it is included everywhere
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
#define AMREX_DEVICE_COMPILE 1
#endif

#else

#define AMREX_GPU_HOST
#define AMREX_GPU_DEVICE
#define AMREX_GPU_GLOBAL
#define AMREX_GPU_HOST_DEVICE
#define AMREX_GPU_DEVICE_MANAGED

#endif

#endif
