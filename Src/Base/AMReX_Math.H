#ifndef AMREX_MATH_H_
#define AMREX_MATH_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <cmath>
#include <cstdlib>

#ifdef AMREX_USE_DPCPP
#include <CL/sycl.hpp>
namespace sycl = cl::sycl;
#endif

namespace amrex { inline namespace disabled {
    // If it is inside namespace amrex, or amrex namespace is imported with using namespace amrex or
    // amrex::disabled, unqualified abs functions are disabled with a compile time error such as,
    // call of overload abs(int&) is ambiguous, or a link time error such as, undefined reference to
    // `amrex::disabled::abs(double)'.  To fix it, one can use `std::abs` or `amrex::Math::abs`.
    // The latter works in both host and device functions, whereas `std::abs` does not currently
    // work on device with HIP and DPC++.
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
//template <typename T> T abs (T a) { return sycl::abs(a); } // The compiler seems to have trouble with this on DG1.
template <typename T> T abs (T a) { return (a >= T(0)) ? a : -a; }

using sycl::ceil;
using sycl::copysign;
using sycl::floor;
using sycl::round;

#else

using std::abs;
using std::ceil;
using std::copysign;
using std::floor;
using std::round;

#endif

}}

#endif
