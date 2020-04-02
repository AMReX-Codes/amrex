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
#else
#define AMREX_GPU_MANAGED __managed__
#define AMREX_GPU_DEVICE_MANAGED __device__ __managed__
#endif

#else

#define AMREX_GPU_HOST
#define AMREX_GPU_DEVICE
#define AMREX_GPU_GLOBAL
#define AMREX_GPU_HOST_DEVICE
#define AMREX_GPU_CONSTANT
#define AMREX_GPU_MANAGED
#define AMREX_GPU_DEVICE_MANAGED

#endif

#define AMREX_DEVICE_COMPILE (__CUDA_ARCH__ || __HIP_DEVICE_COMPILE__ || __SYCL_DEVICE_ONLY__)

#ifdef AMREX_USE_DPCPP
#if 0
// xxxxx DPCPP todo: AMREX_REQUIRE_SUBGROUP_SIZE causes "OpenCL API failed. OpenCL API returns: -46 (CL_INVALID_KERNEL_NAME) -46 (CL_INVALID_KERNEL_NAME)"
# define AMREX_REQUIRE_SUBGROUP_SIZE(x) \
  _Pragma("clang diagnostic push") \
  _Pragma("clang diagnostic ignored \"-Wattributes\"") \
  [[cl::intel_reqd_sub_group_size(x)]] \
  _Pragma("clang diagnostic pop")
#else
#define AMREX_REQUIRE_SUBGROUP_SIZE(x)
#endif
#endif // AMREX_USE_DPCPP

#endif
