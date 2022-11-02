#ifndef AMREX_MATH_H_
#define AMREX_MATH_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <cmath>
#include <cstdlib>
#include <type_traits>
#include <utility>

#ifdef AMREX_USE_DPCPP
#include <CL/sycl.hpp>
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
using sycl::isfinite;
using sycl::isinf;

#else

using std::abs;
using std::ceil;
using std::copysign;
using std::floor;
using std::round;
using std::isfinite;
using std::isinf;

#endif

template <typename T>
constexpr std::enable_if_t<std::is_floating_point<T>::value,T> pi ()
{
    return T(3.1415926535897932384626433832795029L);
}

//! Return cos(x*pi) given x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double cospi (double x)
{
#if defined(AMREX_USE_DPCPP)
    return sycl::cospi(x);
#elif defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    return ::cospi(x);
#else
    return std::cos(pi<double>()*x);
#endif
}

//! Return cos(x*pi) given x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
float cospi (float x)
{
#if defined(AMREX_USE_DPCPP)
    return sycl::cospi(x);
#elif defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    return ::cospif(x);
#else
    return std::cos(pi<float>()*x);
#endif
}

//! Return sin(x*pi) given x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double sinpi (double x)
{
#if defined(AMREX_USE_DPCPP)
    return sycl::sinpi(x);
#elif defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    return ::sinpi(x);
#else
    return std::sin(pi<double>()*x);
#endif
}

//! Return sin(x*pi) given x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
float sinpi (float x)
{
#if defined(AMREX_USE_DPCPP)
    return sycl::sinpi(x);
#elif defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    return ::sinpif(x);
#else
    return std::sin(pi<float>()*x);
#endif
}

//! Return sine and cosine of given number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::pair<double,double> sincos (double x)
{
    std::pair<double,double> r;
#if defined(AMREX_USE_DPCPP)
    r.first = sycl::sincos(x, &r.second);
#elif defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(_GNU_SOURCE)
    ::sincos(x, &r.first, &r.second);
#else
    r.first  = std::sin(x);
    r.second = std::cos(x);
#endif
    return r;
}

//! Return sine and cosine of given number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::pair<float,float> sincos (float x)
{
    std::pair<float,float> r;
#if defined(AMREX_USE_DPCPP)
    r.first = sycl::sincos(x, &r.second);
#elif defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(_GNU_SOURCE)
    ::sincosf(x, &r.first, &r.second);
#else
    r.first  = std::sin(x);
    r.second = std::cos(x);
#endif
    return r;
}

//! Return sin(pi*x) and cos(pi*x) given x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::pair<double,double> sincospi (double x)
{
    std::pair<double,double> r;
#if defined(AMREX_USE_DPCPP)
    r.first = sycl::sincos(pi<double>()*x, &r.second);
#elif defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    ::sincospi(x, &r.first, &r.second);
#else
    r.first  = std::sin(pi<double>()*x);
    r.second = std::cos(pi<double>()*x);
#endif
    return r;
}

//! Return sin(pi*x) and cos(pi*x) given x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::pair<float,float> sincospi (float x)
{
    std::pair<float,float> r;
#if defined(AMREX_USE_DPCPP)
    r.first = sycl::sincos(pi<float>()*x, &r.second);
#elif defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    ::sincospif(x, &r.first, &r.second);
#else
    r.first  = std::sin(pi<float>()*x);
    r.second = std::cos(pi<float>()*x);
#endif
    return r;
}

}}

#endif
