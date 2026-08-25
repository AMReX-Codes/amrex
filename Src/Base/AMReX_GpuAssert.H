#ifndef AMREX_GPU_ASSERT_H_
#define AMREX_GPU_ASSERT_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_GpuAssert.H
 * \brief Assertions for GPU device code.
 *
 * \def AMREX_DEVICE_ASSERT
 * \brief Assert \p flag in CUDA, HIP, or SYCL device code.
 *
 * The macro maps to \c assert and is defined only for those GPU backends on
 * non-Apple platforms.
 */

#if !defined(__APPLE__)

#include <cassert>

#ifdef AMREX_USE_SYCL
#  include <sycl/sycl.hpp>
#endif

#if defined(AMREX_USE_SYCL)
#  define AMREX_DEVICE_ASSERT(flag) assert(flag)
#elif defined(AMREX_USE_CUDA)
#  define AMREX_DEVICE_ASSERT(flag) assert(flag)
#elif defined(AMREX_USE_HIP)
#  define AMREX_DEVICE_ASSERT(flag) assert(flag)
#endif

#endif  // !defined(__APPLE__)
#endif  // AMREX_GPU_ASSERT_H_
