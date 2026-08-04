#ifndef AMREX_COMPENSATEDSUM_H_
#define AMREX_COMPENSATEDSUM_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

namespace amrex
{

// Kahan compensation relies on (t - sum) - y evaluating to the rounding
// error lost by sum + y, rather than to the exact-arithmetic value of 0.
// Optimizations that assume floating-point addition is associative (e.g.
// -ffast-math, -fassociative-math) or that fuse operations across the
// subtractions (aggressive FP contraction) are free to simplify that
// expression to 0, which silently turns compensatedAdd into a plain,
// uncompensated sum. Most compilers only enable such optimizations when a
// user opts in (e.g. AMReX_FASTMATH, off by default -- see
// Tools/CMake/AMReXOptions.cmake). Since this function's correctness must not
// depend on the flags used by its translation unit, its body is pinned to
// strict floating-point semantics regardless.
//
// This function is deliberately NOT AMREX_FORCE_INLINE. Forcing
// inlining across an optimization-setting boundary like the one below
// causes both compilers to inline the -ffast-math-optimized call-site
// version instead of honoring this function's own pragma-scoped settings,
// silently defeating the protection below.
#if defined(__CUDACC__) || defined(__HIPCC__)
#pragma STDC FP_CONTRACT OFF
#endif

#if defined(__clang__)
#pragma clang fp reassociate(off) contract(off)
#elif defined(__GNUC__)
#pragma GCC push_options
#pragma GCC optimize ("no-fast-math")
#endif

/**
 * \brief Add a value to a running sum with Kahan compensation.
 *
 * Callers own both the accumulated sum and the compensation storage. Both
 * should be initialized before the first call, usually to zero.
 */
template <typename T, typename U>
AMREX_GPU_HOST_DEVICE
void compensatedAdd (T& sum, T& compensation, U value) noexcept
{
    T const y = static_cast<T>(value) - compensation;
    T const t = sum + y;
    compensation = (t - sum) - y;
    sum = t;
}

#if defined(__clang__)
#pragma clang fp reassociate(on) contract(on)
#elif defined(__GNUC__)
#pragma GCC pop_options
#endif

}

#endif
