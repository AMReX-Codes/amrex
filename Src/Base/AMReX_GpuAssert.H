#ifndef AMREX_GPU_ASSERT_H_
#define AMREX_GPU_ASSERT_H_

#if !defined(__APPLE__)

#include <cassert>

#ifdef AMREX_USE_DPCPP
#include <CL/sycl.hpp>
#endif

#if defined(AMREX_USE_DPCPP)
// assert(0) is supposed to work. But it does not yet.  abort() works
// by generating the following error:
//      amrex::Abort::0::The program was built for 1 devices
//      Build program log for 'Intel(R) Gen9 HD Graphics NEO':
//      error: undefined reference to `abort()'
//      error: backend compiler failed build.
//       0 (CL_SUCCESS) !!!
//#  define AMREX_DEVICE_ASSERT(flag) assert(flag)
#  define AMREX_DEVICE_ASSERT(flag) abort()
#elif defined(AMREX_USE_CUDA)
#  define AMREX_DEVICE_ASSERT(flag) assert(flag)
#elif defined(AMREX_USE_HIP)
// HIP FIX HERE - assert
#  define AMREX_DEVICE_ASSERT(flag)
#endif

#endif  // !defined(__APPLE__)
#endif  // AMREX_GPU_ASSERT_H_
