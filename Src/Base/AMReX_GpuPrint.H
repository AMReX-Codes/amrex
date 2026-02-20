#ifndef AMREX_GPU_PRINT_H_
#define AMREX_GPU_PRINT_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_GpuPrint.H
 * \brief \c AMREX_DEVICE_PRINTF macro for printing from GPU device code.
 *
 * \def AMREX_DEVICE_PRINTF
 * \brief \c printf-style formatted print callable from GPU device code.
 *
 * Maps to \c printf (CPU fallback), \c sycl::ext::oneapi::experimental::printf (SYCL),
 * or the native CUDA/HIP \c printf depending on the build configuration.
 */

#include <cstdio>

#ifdef AMREX_USE_SYCL
#  include <sycl/sycl.hpp>
#endif

#if defined(AMREX_USE_SYCL)
#if defined(__SYCL_DEVICE_ONLY__)
#  define AMREX_DEVICE_PRINTF(format,...) { \
      static const __attribute__((opencl_constant)) char amrex_i_format[] = format ; \
      sycl::ext::oneapi::experimental::printf(amrex_i_format, __VA_ARGS__); }
#else
#  define AMREX_DEVICE_PRINTF(format,...) { \
      std::printf(format, __VA_ARGS__); }
#endif
#elif defined(AMREX_USE_CUDA)
#  define AMREX_DEVICE_PRINTF(...) std::printf(__VA_ARGS__);
#elif defined(AMREX_USE_HIP)
#  define AMREX_DEVICE_PRINTF(...) ::printf(__VA_ARGS__);
#else
#  define AMREX_DEVICE_PRINTF(format,...) { \
      std::printf(format, __VA_ARGS__); }
#endif

#endif  // AMREX_GPU_PRINT_H_
