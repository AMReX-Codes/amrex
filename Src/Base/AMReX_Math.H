#ifndef AMREX_MATH_H_
#define AMREX_MATH_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <cmath>
#include <cstdlib>
#include <type_traits>
#include <utility>

#ifdef AMREX_USE_SYCL
#  include <sycl/sycl.hpp>
#endif

namespace amrex { inline namespace disabled {
    // If it is inside namespace amrex, or amrex namespace is imported with using namespace amrex or
    // amrex::disabled, unqualified abs functions are disabled with a compile time error such as,
    // call of overload abs(int&) is ambiguous, or a link time error such as, undefined reference to
    // `amrex::disabled::abs(double)'.  To fix it, one can use `std::abs` or `amrex::Math::abs`.
    // The latter works in both host and device functions, whereas `std::abs` does not currently
    // work on device with HIP and SYCL.
    AMREX_GPU_HOST_DEVICE double abs (double);
    AMREX_GPU_HOST_DEVICE float abs (float);
    AMREX_GPU_HOST_DEVICE long double abs (long double);
    AMREX_GPU_HOST_DEVICE int abs (int);
    AMREX_GPU_HOST_DEVICE long abs (long);
    AMREX_GPU_HOST_DEVICE long long abs (long long);
}}

namespace amrex::Math {

// Since Intel's SYCL compiler now supports the following std functions on device,
// one no longer needs to use amrex::Math::abs, etc.  They are kept here for
// backward compatibility.

using std::abs;
using std::ceil;
using std::copysign;
using std::floor;
using std::round;

// However, since Intel's SYCL compiler is very aggressive with fast floating
// point optimisations, the following must be kept, as using the std functions
// always evaluates to false (even at -O1).

#ifdef AMREX_USE_SYCL

using sycl::isfinite;
using sycl::isinf;

#else

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
#if defined(AMREX_USE_SYCL)
    return sycl::cospi(x);
#else
    AMREX_IF_ON_DEVICE(( return ::cospi(x); ))
    AMREX_IF_ON_HOST(( return std::cos(pi<double>()*x); ))
#endif
}

//! Return cos(x*pi) given x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
float cospi (float x)
{
#if defined(AMREX_USE_SYCL)
    return sycl::cospi(x);
#else
    AMREX_IF_ON_DEVICE(( return ::cospif(x); ))
    AMREX_IF_ON_HOST(( return std::cos(pi<float>()*x); ))
#endif
}

//! Return sin(x*pi) given x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double sinpi (double x)
{
#if defined(AMREX_USE_SYCL)
    return sycl::sinpi(x);
#else
    AMREX_IF_ON_DEVICE(( return ::sinpi(x); ))
    AMREX_IF_ON_HOST(( return std::sin(pi<double>()*x); ))
#endif
}

//! Return sin(x*pi) given x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
float sinpi (float x)
{
#if defined(AMREX_USE_SYCL)
    return sycl::sinpi(x);
#else
    AMREX_IF_ON_DEVICE(( return ::sinpif(x); ))
    AMREX_IF_ON_HOST(( return std::sin(pi<float>()*x); ))
#endif
}

namespace detail {
    AMREX_FORCE_INLINE void sincos (double x, double* sinx, double* cosx) {
#if defined(_GNU_SOURCE) && !defined(__APPLE__)
        ::sincos(x, sinx, cosx);
#else
        *sinx = std::sin(x);
        *cosx = std::cos(x);
#endif
    }

    AMREX_FORCE_INLINE void sincosf (float x, float* sinx, float* cosx) {
#if defined(_GNU_SOURCE) && !defined(__APPLE__)
        ::sincosf(x, sinx, cosx);
#else
        *sinx = std::sin(x);
        *cosx = std::cos(x);
#endif
    }
}

//! Return sine and cosine of given number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::pair<double,double> sincos (double x)
{
    std::pair<double,double> r;
#if defined(AMREX_USE_SYCL)
    r.first = sycl::sincos(x, sycl::private_ptr<double>(&r.second));
#else
    AMREX_IF_ON_DEVICE(( ::sincos(x, &r.first, &r.second); ))
    AMREX_IF_ON_HOST(( detail::sincos(x, &r.first, &r.second); ))
#endif
    return r;
}

//! Return sine and cosine of given number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::pair<float,float> sincos (float x)
{
    std::pair<float,float> r;
#if defined(AMREX_USE_SYCL)
    r.first = sycl::sincos(x, sycl::private_ptr<float>(&r.second));
#else
    AMREX_IF_ON_DEVICE(( ::sincosf(x, &r.first, &r.second); ))
    AMREX_IF_ON_HOST(( detail::sincosf(x, &r.first, &r.second); ))
#endif
    return r;
}

//! Return sin(pi*x) and cos(pi*x) given x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::pair<double,double> sincospi (double x)
{
    std::pair<double,double> r;
#if defined(AMREX_USE_SYCL)
    r = sincos(pi<double>()*x);
#else
    AMREX_IF_ON_DEVICE(( ::sincospi(x, &r.first, &r.second); ))
    AMREX_IF_ON_HOST(( r = sincos(pi<double>()*x); ))
#endif
    return r;
}

//! Return sin(pi*x) and cos(pi*x) given x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::pair<float,float> sincospi (float x)
{
    std::pair<float,float> r;
#if defined(AMREX_USE_SYCL)
    r = sincos(pi<float>()*x);
#else
    AMREX_IF_ON_DEVICE(( ::sincospif(x, &r.first, &r.second); ))
    AMREX_IF_ON_HOST(( r = sincos(pi<float>()*x); ))
#endif
    return r;
}

//! Return pow(x, Power), where Power is an integer known at compile time
template <int Power, typename T,
    typename = typename std::enable_if<!std::is_integral<T>() || Power>=0>::type>
AMREX_FORCE_INLINE
constexpr T powi (T x) noexcept
{
    if constexpr (Power < 0) {
        return T(1)/powi<-Power>(x);
    } else if constexpr (Power == 0) {
        //note: 0^0 is implementation-defined, but most compilers return 1
        return T(1);
    } else if constexpr (Power == 1) {
        return x;
    } else if constexpr (Power == 2) {
        return x*x;
    } else if constexpr (Power%2 == 0) {
        return powi<2>(powi<Power/2>(x));
    } else {
        return x*powi<Power-1>(x);
    }
}

}

#endif
