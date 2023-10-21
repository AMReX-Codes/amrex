#ifndef AMREX_GPU_QUALIFIERS_H_
#define AMREX_GPU_QUALIFIERS_H_
#include <AMReX_Config.H>

#if defined(AMREX_USE_GPU) && !defined(AMREX_USE_SYCL)

#if defined(AMREX_USE_HIP)
#include <hip/hip_runtime.h>
#endif

#if defined(AMREX_USE_CUDA) && (defined(AMREX_CXX_PGI) || defined(AMREX_CXX_NVHPC))
#include <nv/target>
#define AMREX_IF_ON_DEVICE(CODE) NV_IF_TARGET(NV_IS_DEVICE, CODE)
#define AMREX_IF_ON_HOST(CODE) NV_IF_TARGET(NV_IS_HOST, CODE)
#endif

#define AMREX_GPU_HOST __host__
#define AMREX_GPU_DEVICE __device__
#define AMREX_GPU_GLOBAL __global__
#define AMREX_GPU_HOST_DEVICE __host__ __device__
#define AMREX_GPU_CONSTANT __constant__

#define AMREX_GPU_MANAGED __managed__
#define AMREX_GPU_DEVICE_MANAGED __device__ __managed__

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

// Remove surrounding parentheses if present
#define AMREX_IMPL_STRIP_PARENS(X) AMREX_IMPL_ESC(AMREX_IMPL_ISH X)
#define AMREX_IMPL_ISH(...) AMREX_IMPL_ISH __VA_ARGS__
#define AMREX_IMPL_ESC(...) AMREX_IMPL_ESC_(__VA_ARGS__)
#define AMREX_IMPL_ESC_(...) AMREX_IMPL_VAN_##__VA_ARGS__
#define AMREX_IMPL_VAN_AMREX_IMPL_ISH

#if !defined(AMREX_IF_ON_DEVICE) && !defined(AMREX_IF_ON_HOST)
#if (defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)) ||         \
    (defined(AMREX_USE_HIP) && defined(__HIP_DEVICE_COMPILE__)) || \
    (defined(AMREX_USE_SYCL) && defined(__SYCL_DEVICE_ONLY__))
#define AMREX_IF_ON_DEVICE(CODE) \
  { AMREX_IMPL_STRIP_PARENS(CODE) }
#define AMREX_IF_ON_HOST(CODE) \
  {}
#else
#define AMREX_IF_ON_DEVICE(CODE) \
  {}
#define AMREX_IF_ON_HOST(CODE) \
  { AMREX_IMPL_STRIP_PARENS(CODE) }
#endif
#endif

#ifdef AMREX_USE_SYCL
# include <sycl/sycl.hpp>
#endif

#endif
