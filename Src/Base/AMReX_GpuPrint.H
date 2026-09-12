#ifndef AMREX_GPU_PRINT_H_
#define AMREX_GPU_PRINT_H_
#include <AMReX_Config.H>

#include <cstdio>

#ifdef AMREX_USE_SYCL
#  include <sycl/sycl.hpp>
#endif

#if defined(AMREX_USE_SYCL)
#  define AMREX_DEVICE_PRINTF(...) \
      sycl::ext::oneapi::experimental::printf(__VA_ARGS__);
#elif defined(AMREX_USE_CUDA)
#  define AMREX_DEVICE_PRINTF(...) std::printf(__VA_ARGS__);
#elif defined(AMREX_USE_HIP)
#  define AMREX_DEVICE_PRINTF(...) ::printf(__VA_ARGS__);
#else
#  define AMREX_DEVICE_PRINTF(...) std::printf(__VA_ARGS__);
#endif

#endif  // AMREX_GPU_PRINT_H_
