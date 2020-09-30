#ifndef AMREX_GPU_ASSERT_H_
#define AMREX_GPU_ASSERT_H_

#if !defined(__APPLE__)

#include <cassert>

#ifdef AMREX_USE_DPCPP
#include <CL/sycl.hpp>
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
