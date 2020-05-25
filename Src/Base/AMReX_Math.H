#ifndef AMREX_MATH_H_
#define AMREX_MATH_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <cmath>

#ifdef AMREX_USE_DPCPP
#include <CL/sycl.hpp>
namespace sycl = cl::sycl;
#endif

namespace amrex { inline namespace disabled {
    // disable unqualified abs function
    AMREX_GPU_HOST_DEVICE double abs (double);
    AMREX_GPU_HOST_DEVICE float abs (float);
    AMREX_GPU_HOST_DEVICE long double abs (long double);
    AMREX_GPU_HOST_DEVICE int abs (int);
    AMREX_GPU_HOST_DEVICE long abs (long);
    AMREX_GPU_HOST_DEVICE long long abs (long long);
}}

namespace amrex { namespace Math {

#ifdef AMREX_USE_DPCPP

//List of unsupported math function can be found at
// https://github.com/intel/llvm/blob/sycl/sycl/doc/extensions/C-CXX-StandardLibrary/C-CXX-StandardLibrary.rst

//using sycl::abs;
// We have to do this because somehow sycl::abs(int) return unsigned int.
template <typename T> T abs (T a) { return sycl::abs(a); }

using sycl::ceil;
using sycl::copysign;
using sycl::floor;
using sycl::round;

using std::sqrt;

#elif defined (AMREX_USE_HIP)

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double abs (double a) { return fabs(a); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
float abs (float a) { return fabsf(a); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int abs (int a) { return ::abs(a); }

using ::ceil;
using ::copysign;
using ::floor;
using ::round;

using ::sqrt;

#else

using std::abs;
using std::ceil;
using std::copysign;
using std::floor;
using std::round;

using std::sqrt;

#endif

}}

#endif
