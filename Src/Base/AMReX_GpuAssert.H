#ifndef AMREX_GPU_ASSERT_H_
#define AMREX_GPU_ASSERT_H_
#include <AMReX_Config.H>

#if !defined(__APPLE__)

#include <cassert>

#ifdef AMREX_USE_DPCPP
#  if defined(__INTEL_LLVM_COMPILER) && (__INTEL_LLVM_COMPILER < 20230100)
#    include <CL/sycl.hpp>
#  else
#    include <sycl/sycl.hpp>
#  endif
#endif

#if defined(AMREX_USE_DPCPP)
#  define AMREX_DEVICE_ASSERT(flag) assert(flag)
#elif defined(AMREX_USE_CUDA)
#  define AMREX_DEVICE_ASSERT(flag) assert(flag)
#elif defined(AMREX_USE_HIP)
#  define AMREX_DEVICE_ASSERT(flag) assert(flag)
#endif

#endif  // !defined(__APPLE__)
#endif  // AMREX_GPU_ASSERT_H_
