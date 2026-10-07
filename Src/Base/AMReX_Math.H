#ifndef AMREX_MATH_H_
#define AMREX_MATH_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <AMReX_INT.H>
#include <AMReX_SIMD.H>
#include <AMReX_REAL.H>
#include <bit>
#include <concepts>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <numbers>
#include <type_traits>
#include <utility>

#ifdef AMREX_USE_SYCL
#  include <sycl/sycl.hpp>
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

// Under fast math (e.g., -ffast-math, -ffinite-math-only, /fp:fast, and Intel's
// SYCL compiler by default), std::isnan etc. may always evaluate to false.  The
// functions below work in all modes.  In fast-math modes, they test the bit
// pattern, hidden from the optimizer so that it cannot fold the test away.

/// \cond DOXYGEN_IGNORE
// Clang: no macro tells if -fno-honor-nans or -fno-honor-infinities is on.
// _WIN32: std::isnan is host only for nvcc with MSVC.
#if (defined(__FINITE_MATH_ONLY__) && (__FINITE_MATH_ONLY__ == 1)) \
    || defined(__FAST_MATH__) || defined(_M_FP_FAST) || defined(_WIN32) \
    || defined(__clang__)
#  define AMREX_MATH_ISNAN_BITS 1
#endif

namespace detail {
    //! Return v, opaque to the optimizer under fast math.  Otherwise, clang
    //! assumes that the result of a fast-math operation is not NaN or inf, and
    //! folds a bit test of it to a constant.
    template <typename U>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    U fp_opaque (U v) noexcept
    {
#if defined(AMREX_MATH_ISNAN_BITS) && (defined(__GNUC__) || defined(__clang__))
#  if defined(__HIP_DEVICE_COMPILE__)
        __asm__("" : "+v"(v));
#  else
        AMREX_IF_ON_HOST((__asm__("" : "+r"(v));))
#  endif
#endif
        return v;
    }

    template <typename T>
    using fp_uint_t = std::conditional_t<sizeof(T) == 8, std::uint64_t, std::uint32_t>;

    //! Exponent bits of T, which are also the bits of +inf.
    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr fp_uint_t<T> fp_exp_mask () noexcept
    {
        return (sizeof(T) == 8) ? fp_uint_t<T>(0x7ff0'0000'0000'0000ULL)
                                : fp_uint_t<T>(0x7f80'0000U);
    }

    //! |x| as bits: NaN if > fp_exp_mask(), inf if ==, finite if <.
    template <bool opaque = true, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    fp_uint_t<T> fp_abs_bits (T x) noexcept
    {
        using U = fp_uint_t<T>;
        auto bits = std::bit_cast<U>(x);
        if constexpr (opaque) { bits = fp_opaque(bits); }
        return bits & (~U(0) >> 1);
    }

    template <typename T>
    concept FloatOrDouble = std::same_as<T,float> || std::same_as<T,double>;
}
/// \endcond

//! Return true if x is NaN.  Works under fast math.
template <typename T>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool isnan (T x) noexcept
{
#if defined(AMREX_USE_SYCL)
    return sycl::isnan(x);
#else
#if defined(AMREX_MATH_ISNAN_BITS)
    if constexpr (detail::FloatOrDouble<T>) {
        return detail::fp_abs_bits(x) > detail::fp_exp_mask<T>();
    } else
#endif
    {
        return std::isnan(x);
    }
#endif
}

//! Return true if x is +/-infinity.  Works under fast math.
template <typename T>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool isinf (T x) noexcept
{
#if defined(AMREX_USE_SYCL)
    return sycl::isinf(x);
#else
#if defined(AMREX_MATH_ISNAN_BITS)
    if constexpr (detail::FloatOrDouble<T>) {
        return detail::fp_abs_bits(x) == detail::fp_exp_mask<T>();
    } else
#endif
    {
        return std::isinf(x);
    }
#endif
}

//! Return true if x is neither NaN nor infinity.  Works under fast math.
template <typename T>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool isfinite (T x) noexcept
{
#if defined(AMREX_USE_SYCL)
    return sycl::isfinite(x);
#else
#if defined(AMREX_MATH_ISNAN_BITS)
    if constexpr (detail::FloatOrDouble<T>) {
        return detail::fp_abs_bits(x) < detail::fp_exp_mask<T>();
    } else
#endif
    {
        return std::isfinite(x);
    }
#endif
}

/// \cond DOXYGEN_IGNORE
namespace detail {
    //! |x|, with NaN mapped to +inf so that a max reduction cannot drop it.
    //! Create it before the loop.  It holds an inf hidden from the optimizer,
    //! because under fast math clang assumes that a result is never inf.
    template <typename T>
    struct AbsNanToInf
    {
        T inf = make_inf();

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T operator() (T x) const noexcept
        {
            if constexpr (FloatOrDouble<T>) {
                auto const a = fp_abs_bits<false>(x);
                return (a < fp_exp_mask<T>()) ? std::bit_cast<T>(a) : inf;
            } else if constexpr (std::is_floating_point_v<T>) {
                return amrex::Math::isnan(x) ? inf : std::abs(x);
            } else {
                return std::abs(x);
            }
        }

    private:
        AMREX_GPU_HOST_DEVICE
        static T make_inf () noexcept
        {
            if constexpr (FloatOrDouble<T>) {
                return std::bit_cast<T>(fp_opaque(fp_exp_mask<T>()));
            } else if constexpr (std::is_floating_point_v<T>) {
                return std::numeric_limits<T>::infinity();
            } else {
                return T(0);
            }
        }
    };
}
/// \endcond

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
}
/// \endcond

#ifdef AMREX_USE_SIMD
//! Return sine and cosine of given number
template<typename T_Real>
requires (amrex::simd::stdx::is_simd_v<T_Real>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::pair<T_Real,T_Real> sincos (T_Real x)
{
    using namespace amrex::simd::stdx;
    std::pair<T_Real,T_Real> r;
    r.first = sin(x);
    r.second = cos(x);
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
std::pair<T_Real,T_Real> sincospi (T_Real x)
{
    using namespace amrex::simd::stdx;
    T_Real const px = pi<typename T_Real::value_type>() * x;
    std::pair<T_Real,T_Real> r;
    r.first = sin(px);
    r.second = cos(px);
    return r;
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
    if (std::abs(k) == T(1)) {
        return std::numeric_limits<T>::infinity();
    }
    // Computing K based on DLMF
    // https://dlmf.nist.gov/19.8
    T tol = std::numeric_limits<T>::epsilon();

    T a0 = 1.0;
    T g0 = static_cast<T>(std::sqrt(1.0 - k*k));
    T a = a0;
    T g = g0;

    // Find Arithmetic Geometric mean
    while(std::abs(a0 - g0) > tol) {
        a = static_cast<T>(0.5*(a0 + g0));
        g = std::sqrt(a0 * g0);

        a0 = a;
        g0 = g;
    }

    return static_cast<T>(0.5*pi<T>()/a);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T comp_ellint_2 (T k)
{
    if (std::abs(k) == T(1)) {
        return T(1);
    }
    // Computing E based on DLMF
    // https://dlmf.nist.gov/19.8
    T Kcomp = amrex::Math::comp_ellint_1<T>(k);
    T tol = std::numeric_limits<T>::epsilon();

    // Step Zero
    T a0 = 1.0;
    T g0 = static_cast<T>(std::sqrt(1.0 - k*k));
    T cn = std::sqrt(a0*a0 - g0*g0);

    // Step 1
    int n = 1;
    T a = static_cast<T>(0.5 * (a0 + g0));
    T g = std::sqrt(a0*g0);
    cn = static_cast<T>(0.25*cn*cn/a);

    T sum_val = a*a;
    a0 = a;
    g0 = g;

    while(std::abs(cn*cn) > tol) {
        // Compute coefficients for this iteration
        a = static_cast<T>(0.5 * (a0 + g0));
        g = std::sqrt(a0*g0);
        cn = static_cast<T>(0.25*cn*cn/a);

        n++;
        sum_val = static_cast<T>(sum_val - std::pow(2,n-1)*cn*cn);

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

namespace amrex {
    using Math::isnan;
    using Math::isinf;
    using Math::isfinite;
}

#endif
