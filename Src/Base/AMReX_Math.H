#ifndef AMREX_MATH_H_
#define AMREX_MATH_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <AMReX_INT.H>
#include <AMReX_SIMD.H>
#include <AMReX_REAL.H>
#include <concepts>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <numbers>
#include <type_traits>
#include <utility>

#ifdef AMREX_USE_SYCL
#  include <sycl/sycl.hpp>
#endif

/*
 * Vector math library support for the SIMD math functions in amrex::Math
 *
 * With AMReX_SIMD=ON, the amrex::Math transcendentals below also accept SIMD
 * types. They are implemented as a short loop over the SIMD lanes with a
 * compile-time constant trip count (see amrex::Math::detail::map_lanes).
 * Compilers that know a vector variant of the scalar function called inside
 * that loop collapse the whole loop into a single call into a vector math
 * library, such as glibc's libmvec.
 *
 * Two conditions must hold for that to happen:
 *  - the translation unit must be compiled with -fno-math-errno, so that a math
 *    call is free of side effects and may be vectorized, and
 *  - the vector variants must have been declared. glibc declares them only for
 *    translation units compiled with -ffast-math, which is too broad a hammer
 *    for most codes, so the ones AMReX uses are declared here instead.
 *
 * The CMake option AMReX_SIMD_VECMATH=ON adds the required compiler flags. If
 * either condition is missing, the loop merely unrolls into scalar calls, which
 * is what the SIMD library's own fallback does as well.
 *
 * Note on accuracy: vector math libraries trade accuracy for speed. glibc's
 * libmvec documents a maximum error of 4 ULP, where its scalar routines stay
 * below 1 ULP. Results are also no longer bit-wise identical to a scalar build.
 *
 * @todo Move the SIMD math functions and the declarations below into their own
 *       header, e.g. AMReX_Math_SIMD.H, that AMReX_Math.H includes. Keeping the
 *       vector variant declarations out of a header that is pulled in almost
 *       everywhere makes it clearer where they take effect.
 */
/* AMREX_SIMD_VECMATH is defined when the lane loops below are expected to reach
 * a vector math library. It can also be set from the outside (the CMake option
 * AMReX_SIMD_VECMATH=ON does so for compilers we cannot detect from here, such
 * as clang with -fveclib=libmvec).
 */
#if !defined(AMREX_SIMD_VECMATH) && defined(AMREX_USE_SIMD) && defined(__NO_MATH_ERRNO__) \
    && defined(__GLIBC__) && defined(__x86_64__) && defined(__GNUC__) && !defined(__clang__)
#   include <features.h>
#   ifdef __GLIBC_PREREQ
#       if __GLIBC_PREREQ(2,22)
//! Defined when SIMD transcendentals are expected to lower to vector math library calls
#           define AMREX_SIMD_VECMATH 1
#       endif
#   endif
#endif

#if defined(AMREX_SIMD_VECMATH) && defined(__GLIBC__) && defined(__x86_64__) \
    && defined(__GNUC__) && !defined(__clang__) && !defined(__FAST_MATH__)
#   include <features.h>
#   ifdef __GLIBC_PREREQ
// glibc's own declarations use "notinbranch": no masked variant is provided.
#       define AMREX_VECMATH_FN __attribute__((__simd__("notinbranch")))
#       if __GLIBC_PREREQ(2,22)
extern "C" {
    AMREX_VECMATH_FN double cos   (double) noexcept;
    AMREX_VECMATH_FN double exp   (double) noexcept;
    AMREX_VECMATH_FN double log   (double) noexcept;
    AMREX_VECMATH_FN double pow   (double, double) noexcept;
    AMREX_VECMATH_FN double sin   (double) noexcept;
    AMREX_VECMATH_FN float  cosf  (float) noexcept;
    AMREX_VECMATH_FN float  expf  (float) noexcept;
    AMREX_VECMATH_FN float  logf  (float) noexcept;
    AMREX_VECMATH_FN float  powf  (float, float) noexcept;
    AMREX_VECMATH_FN float  sinf  (float) noexcept;
#           if defined(_GNU_SOURCE) && !defined(__APPLE__)
    AMREX_VECMATH_FN void   sincos  (double, double*, double*) noexcept;
    AMREX_VECMATH_FN void   sincosf (float, float*, float*) noexcept;
#           endif
}
#       endif
#       if __GLIBC_PREREQ(2,35)
extern "C" {
    AMREX_VECMATH_FN double acos   (double) noexcept;
    AMREX_VECMATH_FN double acosh  (double) noexcept;
    AMREX_VECMATH_FN double asin   (double) noexcept;
    AMREX_VECMATH_FN double asinh  (double) noexcept;
    AMREX_VECMATH_FN double atan   (double) noexcept;
    AMREX_VECMATH_FN double atan2  (double, double) noexcept;
    AMREX_VECMATH_FN double atanh  (double) noexcept;
    AMREX_VECMATH_FN double cbrt   (double) noexcept;
    AMREX_VECMATH_FN double cosh   (double) noexcept;
    AMREX_VECMATH_FN double erf    (double) noexcept;
    AMREX_VECMATH_FN double erfc   (double) noexcept;
    AMREX_VECMATH_FN double exp2   (double) noexcept;
    AMREX_VECMATH_FN double expm1  (double) noexcept;
    AMREX_VECMATH_FN double hypot  (double, double) noexcept;
    AMREX_VECMATH_FN double log10  (double) noexcept;
    AMREX_VECMATH_FN double log1p  (double) noexcept;
    AMREX_VECMATH_FN double log2   (double) noexcept;
    AMREX_VECMATH_FN double sinh   (double) noexcept;
    AMREX_VECMATH_FN double tan    (double) noexcept;
    AMREX_VECMATH_FN double tanh   (double) noexcept;
    AMREX_VECMATH_FN float  acosf  (float) noexcept;
    AMREX_VECMATH_FN float  acoshf (float) noexcept;
    AMREX_VECMATH_FN float  asinf  (float) noexcept;
    AMREX_VECMATH_FN float  asinhf (float) noexcept;
    AMREX_VECMATH_FN float  atanf  (float) noexcept;
    AMREX_VECMATH_FN float  atan2f (float, float) noexcept;
    AMREX_VECMATH_FN float  atanhf (float) noexcept;
    AMREX_VECMATH_FN float  cbrtf  (float) noexcept;
    AMREX_VECMATH_FN float  coshf  (float) noexcept;
    AMREX_VECMATH_FN float  erff   (float) noexcept;
    AMREX_VECMATH_FN float  erfcf  (float) noexcept;
    AMREX_VECMATH_FN float  exp2f  (float) noexcept;
    AMREX_VECMATH_FN float  expm1f (float) noexcept;
    AMREX_VECMATH_FN float  hypotf (float, float) noexcept;
    AMREX_VECMATH_FN float  log10f (float) noexcept;
    AMREX_VECMATH_FN float  log1pf (float) noexcept;
    AMREX_VECMATH_FN float  log2f  (float) noexcept;
    AMREX_VECMATH_FN float  sinhf  (float) noexcept;
    AMREX_VECMATH_FN float  tanf   (float) noexcept;
    AMREX_VECMATH_FN float  tanhf  (float) noexcept;
}
#       endif
#       undef AMREX_VECMATH_FN
#   endif
#endif

namespace amrex { // NOLINT(modernize-concat-nested-namespaces)
/// \cond DOXYGEN_IGNORE
inline namespace disabled {
    // If it is inside namespace amrex, or amrex namespace is imported with using namespace amrex or
    // amrex::disabled, unqualified abs functions are disabled with a compile time error such as,
    // call of overload abs(int&) is ambiguous, or a link time error such as, undefined reference to
    // `amrex::disabled::abs(double)'.  To fix it, one can use `std::abs` or `amrex::Math::abs`.
    // We have amrex::Math::abs, because std::abs did not work with HIP and SYCL in the past.
    AMREX_GPU_HOST_DEVICE double abs (double);
    AMREX_GPU_HOST_DEVICE float abs (float);
    AMREX_GPU_HOST_DEVICE long double abs (long double);
    AMREX_GPU_HOST_DEVICE int abs (int);
    AMREX_GPU_HOST_DEVICE long abs (long);
    AMREX_GPU_HOST_DEVICE long long abs (long long);
}
/// \endcond
}

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

template <std::floating_point T>
constexpr T pi ()
{
    return std::numbers::pi_v<T>;
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

/// \cond DOXYGEN_IGNORE
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
#ifdef AMREX_USE_SIMD
    /** Apply a scalar function to every lane of a SIMD variable
     *
     * The lane loop has a compile-time constant trip count and a body that is a
     * single scalar function call. That is the shape an auto-vectorizer can
     * replace with one call into a vector math library (see the note near the
     * top of this file). When it does, the lane buffer never reaches memory;
     * when it does not, the loop unrolls into scalar calls, which is what the
     * SIMD library's own math fallback does as well.
     *
     * @param x the SIMD variable to transform
     * @param f a scalar function applied to one lane
     * @return a SIMD variable holding f(x) for every lane
     */
    template <typename T_Simd, typename F>
    AMREX_FORCE_INLINE
    T_Simd map_lanes (T_Simd const& x, F const& f)
    {
#ifdef AMREX_SIMD_VECMATH
        using T = typename T_Simd::value_type;
        constexpr std::size_t width = T_Simd::size();
        constexpr std::size_t alignment = amrex::simd::stdx::memory_alignment_v<T_Simd>;

        alignas(alignment) T lane[width];
        x.copy_to(lane, amrex::simd::stdx::vector_aligned);

        for (std::size_t i = 0; i < width; ++i) {
            lane[i] = f(lane[i]);
        }

        T_Simd r;
        r.copy_from(lane, amrex::simd::stdx::vector_aligned);
        return r;
#else
        // Without a vector math library to call, the loop above would only add a
        // lane buffer around the same scalar calls. Build the result element by
        // element instead, like the SIMD library's own math fallback does.
        return T_Simd([&] (auto i) { return f(x[i]); });
#endif
    }

    /** Apply a scalar function of two arguments to every lane of two SIMD variables
     *
     * @see map_lanes(T_Simd const&, F const&)
     */
    template <typename T_Simd, typename F>
    AMREX_FORCE_INLINE
    T_Simd map_lanes (T_Simd const& x, T_Simd const& y, F const& f)
    {
#ifdef AMREX_SIMD_VECMATH
        using T = typename T_Simd::value_type;
        constexpr std::size_t width = T_Simd::size();
        constexpr std::size_t alignment = amrex::simd::stdx::memory_alignment_v<T_Simd>;

        alignas(alignment) T lane_x[width];
        alignas(alignment) T lane_y[width];
        x.copy_to(lane_x, amrex::simd::stdx::vector_aligned);
        y.copy_to(lane_y, amrex::simd::stdx::vector_aligned);

        for (std::size_t i = 0; i < width; ++i) {
            lane_x[i] = f(lane_x[i], lane_y[i]);
        }

        T_Simd r;
        r.copy_from(lane_x, amrex::simd::stdx::vector_aligned);
        return r;
#else
        // see map_lanes(T_Simd const&, F const&)
        return T_Simd([&] (auto i) { return f(x[i], y[i]); });
#endif
    }

#endif
}
/// \endcond

#ifdef AMREX_USE_SIMD
//! Return sine and cosine of given number
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::pair<T_Real,T_Real> sincos (T_Real const& x)
{
    // Evaluated as two separate lane loops on purpose: no vector math library
    // offers an auto-vectorizable sincos, because writing both results through
    // pointers stops the vectorizer. Two vector calls still beat one scalar
    // sincos per lane.
    std::pair<T_Real,T_Real> r;
    r.first  = detail::map_lanes(x, [] (auto v) { return std::sin(v); });
    r.second = detail::map_lanes(x, [] (auto v) { return std::cos(v); });
    return r;
}
#endif

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

#ifdef AMREX_USE_SIMD
//! Return sin(pi*x) and cos(pi*x) given x
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::pair<T_Real,T_Real> sincospi (T_Real const& x)
{
    return sincos(T_Real(pi<typename T_Real::value_type>()) * x);
}
#endif

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

/*
 * Transcendental functions
 *
 * Each function comes in a scalar and, with AMReX_SIMD=ON, a SIMD overload, so
 * that a compute kernel written once can be instantiated for both. Call them
 * qualified as amrex::Math::sinh(x): an unqualified call on a SIMD argument
 * resolves to the SIMD library's own (scalar, element-wise) overload through
 * argument-dependent lookup instead.
 *
 * The SIMD overloads may lower to vector math library calls, see the note near
 * the top of this file for the requirements and the accuracy implications.
 */

//! Return the sine of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real sin (T_Real x)
{
    return std::sin(x);
}

#ifdef AMREX_USE_SIMD
//! Return the sine of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real sin (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::sin(v); });
}
#endif

//! Return the cosine of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real cos (T_Real x)
{
    return std::cos(x);
}

#ifdef AMREX_USE_SIMD
//! Return the cosine of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real cos (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::cos(v); });
}
#endif

//! Return the tangent of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real tan (T_Real x)
{
    return std::tan(x);
}

#ifdef AMREX_USE_SIMD
//! Return the tangent of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real tan (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::tan(v); });
}
#endif

//! Return the arc sine of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real asin (T_Real x)
{
    return std::asin(x);
}

#ifdef AMREX_USE_SIMD
//! Return the arc sine of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real asin (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::asin(v); });
}
#endif

//! Return the arc cosine of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real acos (T_Real x)
{
    return std::acos(x);
}

#ifdef AMREX_USE_SIMD
//! Return the arc cosine of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real acos (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::acos(v); });
}
#endif

//! Return the arc tangent of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real atan (T_Real x)
{
    return std::atan(x);
}

#ifdef AMREX_USE_SIMD
//! Return the arc tangent of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real atan (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::atan(v); });
}
#endif

//! Return the hyperbolic sine of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real sinh (T_Real x)
{
    return std::sinh(x);
}

#ifdef AMREX_USE_SIMD
//! Return the hyperbolic sine of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real sinh (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::sinh(v); });
}
#endif

//! Return the hyperbolic cosine of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real cosh (T_Real x)
{
    return std::cosh(x);
}

#ifdef AMREX_USE_SIMD
//! Return the hyperbolic cosine of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real cosh (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::cosh(v); });
}
#endif

//! Return the hyperbolic tangent of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real tanh (T_Real x)
{
    return std::tanh(x);
}

#ifdef AMREX_USE_SIMD
//! Return the hyperbolic tangent of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real tanh (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::tanh(v); });
}
#endif

//! Return the inverse hyperbolic sine of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real asinh (T_Real x)
{
    return std::asinh(x);
}

#ifdef AMREX_USE_SIMD
//! Return the inverse hyperbolic sine of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real asinh (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::asinh(v); });
}
#endif

//! Return the inverse hyperbolic cosine of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real acosh (T_Real x)
{
    return std::acosh(x);
}

#ifdef AMREX_USE_SIMD
//! Return the inverse hyperbolic cosine of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real acosh (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::acosh(v); });
}
#endif

//! Return the inverse hyperbolic tangent of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real atanh (T_Real x)
{
    return std::atanh(x);
}

#ifdef AMREX_USE_SIMD
//! Return the inverse hyperbolic tangent of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real atanh (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::atanh(v); });
}
#endif

//! Return the base-e exponential of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real exp (T_Real x)
{
    return std::exp(x);
}

#ifdef AMREX_USE_SIMD
//! Return the base-e exponential of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real exp (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::exp(v); });
}
#endif

//! Return the base-2 exponential of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real exp2 (T_Real x)
{
    return std::exp2(x);
}

#ifdef AMREX_USE_SIMD
//! Return the base-2 exponential of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real exp2 (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::exp2(v); });
}
#endif

//! Return the base-e exponential minus one of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real expm1 (T_Real x)
{
    return std::expm1(x);
}

#ifdef AMREX_USE_SIMD
//! Return the base-e exponential minus one of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real expm1 (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::expm1(v); });
}
#endif

//! Return the natural logarithm of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real log (T_Real x)
{
    return std::log(x);
}

#ifdef AMREX_USE_SIMD
//! Return the natural logarithm of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real log (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::log(v); });
}
#endif

//! Return the base-2 logarithm of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real log2 (T_Real x)
{
    return std::log2(x);
}

#ifdef AMREX_USE_SIMD
//! Return the base-2 logarithm of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real log2 (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::log2(v); });
}
#endif

//! Return the base-10 logarithm of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real log10 (T_Real x)
{
    return std::log10(x);
}

#ifdef AMREX_USE_SIMD
//! Return the base-10 logarithm of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real log10 (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::log10(v); });
}
#endif

//! Return the natural logarithm of one plus the argument of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real log1p (T_Real x)
{
    return std::log1p(x);
}

#ifdef AMREX_USE_SIMD
//! Return the natural logarithm of one plus the argument of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real log1p (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::log1p(v); });
}
#endif

//! Return the cube root of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real cbrt (T_Real x)
{
    return std::cbrt(x);
}

#ifdef AMREX_USE_SIMD
//! Return the cube root of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real cbrt (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::cbrt(v); });
}
#endif

//! Return the error function of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real erf (T_Real x)
{
    return std::erf(x);
}

#ifdef AMREX_USE_SIMD
//! Return the error function of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real erf (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::erf(v); });
}
#endif

//! Return the complementary error function of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real erfc (T_Real x)
{
    return std::erfc(x);
}

#ifdef AMREX_USE_SIMD
//! Return the complementary error function of every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real erfc (T_Real const& x)
{
    return detail::map_lanes(x, [] (auto v) { return std::erfc(v); });
}
#endif

//! Return x raised to the power y
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real pow (T_Real x, T_Real y)
{
    return std::pow(x, y);
}

#ifdef AMREX_USE_SIMD
//! Return x raised to the power y, for every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real pow (T_Real const& x, T_Real const& y)
{
    return detail::map_lanes(x, y, [] (auto v_x, auto v_y) { return std::pow(v_x, v_y); });
}
#endif

//! Return arc tangent of y/x, using their signs
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real atan2 (T_Real y, T_Real x)
{
    return std::atan2(y, x);
}

#ifdef AMREX_USE_SIMD
//! Return arc tangent of y/x, using their signs, for every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real atan2 (T_Real const& y, T_Real const& x)
{
    return detail::map_lanes(y, x, [] (auto v_y, auto v_x) { return std::atan2(v_y, v_x); });
}
#endif

//! Return square root of x*x + y*y
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real hypot (T_Real x, T_Real y)
{
    return std::hypot(x, y);
}

#ifdef AMREX_USE_SIMD
//! Return square root of x*x + y*y, for every SIMD lane
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real hypot (T_Real const& x, T_Real const& y)
{
    return detail::map_lanes(x, y, [] (auto v_x, auto v_y) { return std::hypot(v_x, v_y); });
}
#endif

//! Return the square root of the given number
template<typename T_Real>
requires (std::is_floating_point_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T_Real sqrt (T_Real x)
{
    return std::sqrt(x);
}

#ifdef AMREX_USE_SIMD
/** Return the square root of every SIMD lane
 *
 * Unlike the transcendentals above, this maps onto a hardware instruction and
 * is handed to the SIMD library directly. No vector math library is involved
 * and the result is correctly rounded.
 */
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real sqrt (T_Real const& x)
{
    return amrex::simd::stdx::sqrt(x);
}

/** Return the absolute value of every SIMD lane
 *
 * @see sqrt for why no vector math library is involved
 */
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_FORCE_INLINE
T_Real abs (T_Real const& x)
{
    return amrex::simd::stdx::abs(x);
}
#endif

//! Return pow(x, Power), where Power is an integer known at compile time
template <int Power, typename T>
requires (!std::integral<T> || Power >= 0)
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

} // namespace amrex::Math

#if defined(AMREX_USE_CUDA)
// Forward-declare libdevice integer-exponent intrinsics. They are part of
// NVIDIA's libdevice and linked at PTX time but are not declared in standard
// CUDA cmath headers outside __CUDACC_RTC__ / _LIBCPP_VERSION blocks.
// See https://docs.nvidia.com/cuda/libdevice-users-guide/__nv_powif.html
extern "C" __device__ float  __nv_powif (float,  int);
extern "C" __device__ double __nv_powi  (double, int);
#endif

namespace amrex::Math {

//! Return pow(x, n) with a runtime integer exponent.
//!
//! On GPU backends this routes to the integer-exponent intrinsic
//! (libdevice __nv_powif / __nv_powi on CUDA, __ocml_pown_f32 /
//! __ocml_pown_f64 on HIP, sycl::pown on SYCL), which is faster than the
//! general pow and avoids the C++11 mixed-arg promotion that would
//! otherwise call the slower double-precision pow(double, double) and
//! narrow on return. On host it casts the int exponent to the
//! floating-point type and calls std::pow.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
float powi (float x, int n) noexcept
{
#if defined(AMREX_USE_SYCL)
    return sycl::pown(x, n);
#elif defined(AMREX_USE_HIP)
    AMREX_IF_ON_DEVICE(( return ::powif(x, n); ))
    AMREX_IF_ON_HOST(( return std::pow(x, static_cast<float>(n)); ))
#elif defined(AMREX_USE_CUDA)
    AMREX_IF_ON_DEVICE(( return __nv_powif(x, n); ))
    AMREX_IF_ON_HOST(( return std::pow(x, static_cast<float>(n)); ))
#else
    return std::pow(x, static_cast<float>(n));
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double powi (double x, int n) noexcept
{
#if defined(AMREX_USE_SYCL)
    return sycl::pown(x, n);
#elif defined(AMREX_USE_HIP)
    AMREX_IF_ON_DEVICE(( return ::powi(x, n); ))
    AMREX_IF_ON_HOST(( return std::pow(x, n); ))
#elif defined(AMREX_USE_CUDA)
    AMREX_IF_ON_DEVICE(( return __nv_powi(x, n); ))
    AMREX_IF_ON_HOST(( return std::pow(x, n); ))
#else
    return std::pow(x, n);
#endif
}

#if defined(AMREX_INT128_SUPPORTED)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::uint64_t umulhi (std::uint64_t a, std::uint64_t b)
{
#if defined(AMREX_USE_SYCL)
    return sycl::mul_hi(a,b);
#else
    AMREX_IF_ON_DEVICE(( return __umul64hi(a, b); ))
    AMREX_IF_ON_HOST((
        auto tmp = amrex::UInt128_t(a) * amrex::UInt128_t(b);
        return std::uint64_t(tmp >> 64);
    ))
#endif
}
#endif

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T comp_ellint_1 (T k)
{
    // Computing K based on DLMF
    // https://dlmf.nist.gov/19.8
    T tol = std::numeric_limits<T>::epsilon();

    T a0 = 1.0;
    T g0 = std::sqrt(1.0 - k*k);
    T a = a0;
    T g = g0;

    // Find Arithmetic Geometric mean
    while(std::abs(a0 - g0) > tol) {
        a = 0.5*(a0 + g0);
        g = std::sqrt(a0 * g0);

        a0 = a;
        g0 = g;
    }

    return 0.5*pi<T>()/a;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T comp_ellint_2 (T k)
{
    // Computing E based on DLMF
    // https://dlmf.nist.gov/19.8
    T Kcomp = amrex::Math::comp_ellint_1<T>(k);
    T tol = std::numeric_limits<T>::epsilon();

    // Step Zero
    T a0 = 1.0;
    T g0 = std::sqrt(1.0 - k*k);
    T cn = std::sqrt(a0*a0 - g0*g0);

    // Step 1
    int n = 1;
    T a = 0.5 * (a0 + g0);
    T g = std::sqrt(a0*g0);
    cn = 0.25*cn*cn/a;

    T sum_val = a*a;
    a0 = a;
    g0 = g;

    while(std::abs(cn*cn) > tol) {
        // Compute coefficients for this iteration
        a = 0.5 * (a0 + g0);
        g = std::sqrt(a0*g0);
        cn = 0.25*cn*cn/a;

        n++;
        sum_val -= std::pow(2,n-1)*cn*cn;

        // Save a and g for next iteration
        a0 = a;
        g0 = g;
    }

    return Kcomp*sum_val;
}

//! Return inverse square root of x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double rsqrt (double x)
{
    double r;
#if defined(AMREX_USE_SYCL)
    AMREX_IF_ON_DEVICE(( r = sycl::rsqrt(x); ))
#else
    AMREX_IF_ON_DEVICE(( r = ::rsqrt(x); ))
#endif
    AMREX_IF_ON_HOST(( r = 1. / std::sqrt(x); ))
    return r;
}

//! Return inverse square root of x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
float rsqrt (float x)
{
    float r;
#if defined(AMREX_USE_SYCL)
    AMREX_IF_ON_DEVICE(( r = sycl::rsqrt(x); ))
#else
    AMREX_IF_ON_DEVICE(( r = ::rsqrtf(x); ))
#endif
    AMREX_IF_ON_HOST(( r = 1.F / std::sqrt(x); ))
    return r;
}

//! Return 10**x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double exp10 (double x)
{
    double r;
#if defined(AMREX_USE_SYCL)
    AMREX_IF_ON_DEVICE(( r = sycl::exp10(x); ))
#else
    AMREX_IF_ON_DEVICE(( r = ::exp10(x); ))
#endif
    AMREX_IF_ON_HOST(( r = std::pow(10.0, x); ))
    return r;
}

//! Return 10**x
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
float exp10 (float x)
{
    float r;
#if defined(AMREX_USE_SYCL)
    AMREX_IF_ON_DEVICE(( r = sycl::exp10(x); ))
#else
    AMREX_IF_ON_DEVICE(( r = ::exp10f(x); ))
#endif
    AMREX_IF_ON_HOST(( r = std::pow(10.0F, x); ))
    return r;
}

/***************************************************************************************************
 * Copyright (c) 2017 - 2024 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **************************************************************************************************/

/////////////////////////////////////////////////////////////////////////////////////////////////

/// Object to encapsulate the fast division+modulus operation for 64b integer division.
///
/// Example:
///
///
///   uint64_t quotient, remainder, dividend, divisor;
///
///   FastDivmodU64 divmod(divisor);
///
///   divmod(quotient, remainder, dividend);
///
///   // quotient = (dividend / divisor)
///   // remainder = (dividend % divisor)
///
struct FastDivmodU64
{
    std::uint64_t divisor = 0;

#ifdef AMREX_INT128_SUPPORTED
    std::uint64_t multiplier = 1U;
    unsigned int shift_right = 0;
    unsigned int round_up = 0;

    //
    // Static methods
    //

    /// Computes b, where 2^b is the greatest power of two that is less than or equal to x
    static std::uint32_t integer_log2 (std::uint64_t x)
    {
        std::uint32_t n = 0;
        while (x >>= 1) {
            ++n;
        }
        return n;
    }

    /// Construct the FastDivmod object, in host code only
    ///
    /// This precomputes some values based on the divisor and is computationally expensive.
    FastDivmodU64 (std::uint64_t divisor_)
        : divisor(divisor_)
    {
        if (divisor) {
            shift_right = integer_log2(divisor);

            if ((divisor & (divisor - 1)) == 0) {
                multiplier = 0;
            }
            else {
                std::uint64_t power_of_two = (std::uint64_t(1) << shift_right);
                auto n = amrex::UInt128_t(power_of_two) << 64;
                std::uint64_t multiplier_lo = n / divisor;
                n += power_of_two;
                multiplier = n / divisor;
                round_up = (multiplier_lo == multiplier ? 1 : 0);
            }
        }
    }

#else

    FastDivmodU64 (std::uint64_t divisor_) : divisor(divisor_) {}

#endif

    /// Default construct an invalid FastDivmodU64
    FastDivmodU64 () = default;

    /// Returns the quotient of floor(dividend / divisor)
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    std::uint64_t divide (std::uint64_t dividend) const
    {
#if defined(AMREX_INT128_SUPPORTED)
        auto x = dividend;
        if (multiplier) {
            x = amrex::Math::umulhi(dividend + round_up, multiplier);
        }
        return (x >> shift_right);
#else
        return dividend / divisor;
#endif
    }

    /// Computes the remainder given a computed quotient and dividend
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    std::uint64_t modulus (std::uint64_t quotient, std::uint64_t dividend) const
    {
        return dividend - quotient * divisor;
    }

    /// Returns the quotient of floor(dividend / divisor) and computes the remainder
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    std::uint64_t divmod (std::uint64_t &remainder, std::uint64_t dividend) const
    {
        auto quotient = divide(dividend);
        remainder = modulus(quotient, dividend);
        return quotient;
    }

    /// Computes integer division and modulus using precomputed values. This is computationally
    /// inexpensive.
    AMREX_GPU_HOST_DEVICE
    void operator() (std::uint64_t &quotient, std::uint64_t &remainder, std::uint64_t dividend) const
    {
        quotient = divmod(remainder, dividend);
    }
};

}

#endif
