#ifndef AMREX_GPU_PRINT_H_
#define AMREX_GPU_PRINT_H_

#if !defined(__APPLE__)

#include <cstdio>

#ifdef AMREX_USE_DPCPP
#include <CL/sycl.hpp>
#endif

#if defined(AMREX_USE_DPCPP)
#  define AMREX_DEVICE_PRINTF(format,...) { \
      static const __attribute__((opencl_constant)) char amrex_i_format[] = format ; \
      cl::sycl::intel::experimental::printf(amrex_i_format, __VA_ARGS__); }
#elif defined(AMREX_USE_CUDA)
#  define AMREX_DEVICE_PRINTF(...) std::printf(__VA_ARGS__);
#elif defined(AMREX_USE_HIP)
#  define AMREX_DEVICE_PRINTF(...) ::printf(__VA_ARGS__);
#endif

#endif  // !defined(__APPLE__)
#endif  // AMREX_GPU_PRINT_H_

