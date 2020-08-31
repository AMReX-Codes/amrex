#ifndef AMREX_GPU_QUALIFIERS_H_
#define AMREX_GPU_QUALIFIERS_H_

#if defined(AMREX_USE_GPU) && !defined(AMREX_USE_DPCPP)

#if defined(AMREX_USE_HIP)
#include <hip/hip_runtime.h>
#endif

#define AMREX_GPU_HOST __host__
#define AMREX_GPU_DEVICE __device__
#define AMREX_GPU_GLOBAL __global__
#define AMREX_GPU_HOST_DEVICE __host__ __device__
#define AMREX_GPU_CONSTANT __constant__

#if defined(AMREX_USE_HIP) && !defined(AMREX_USE_CUDA)
#define AMREX_GPU_MANAGED __device__
#define AMREX_GPU_DEVICE_MANAGED __device__
#define AMREX_GPU_EXTERNAL __attribute__((weak))
#else
#define AMREX_GPU_MANAGED __managed__
#define AMREX_GPU_DEVICE_MANAGED __device__ __managed__
#define AMREX_GPU_EXTERNAL
#endif

#else

#define AMREX_GPU_HOST
#define AMREX_GPU_DEVICE
#define AMREX_GPU_GLOBAL
#define AMREX_GPU_HOST_DEVICE
#define AMREX_GPU_CONSTANT
#define AMREX_GPU_MANAGED
#define AMREX_GPU_DEVICE_MANAGED
#ifdef AMREX_USE_DPCPP
#define AMREX_GPU_EXTERNAL SYCL_EXTERNAL
#else
#define AMREX_GPU_EXTERNAL
#endif

#endif

#define AMREX_DEVICE_COMPILE (__CUDA_ARCH__ || __HIP_DEVICE_COMPILE__ || __SYCL_DEVICE_ONLY__)

#ifdef AMREX_USE_DPCPP
# define AMREX_REQUIRE_SUBGROUP_SIZE(x) \
  _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Wattributes\"") \
  [[cl::intel_reqd_sub_group_size(x)]] \
  _Pragma("clang diagnostic pop")
#endif // AMREX_USE_DPCPP

#endif
