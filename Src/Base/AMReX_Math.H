#ifndef AMREX_MATH_H_
#define AMREX_MATH_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <AMReX_INT.H>
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
constexpr std::enable_if_t<std::is_floating_point_v<T>,T> pi ()
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
    typename = std::enable_if_t<!std::is_integral<T>() || Power>=0>>
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
