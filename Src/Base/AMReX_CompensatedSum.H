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
//
// The pragmas below were verified empirically (Compiler Explorer) for GCC,
// Clang, Intel classic icc, Intel icx (clang-based, covered by the clang
// branch), MSVC, clang-cl (defines both _MSC_VER and __clang__; the clang
// branch takes priority and is what actually applies), NVCC host and
// device code, and AMD's HIP/ROCm clang host and device code (also covered
// by the clang branch). NVCC's device compiler was confirmed not to need
// protection: neither --use_fast_math nor --fmad=true collapses this
// expression, unlike host-side -ffast-math. The nv_diag_suppress pair
// exists only to silence NVCC's front end complaining about the (correct,
// but nvcc-unrecognized) GCC/MSVC/Intel pragmas when it parses this header
// while selecting a host-side branch; it has no effect on codegen.
#if defined(__CUDACC__)
#pragma nv_diag_suppress 1675
#endif

#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
#pragma float_control(precise, on, push)
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize ("no-fast-math")
#endif

/**
 * \brief Add a value to a running sum with Kahan compensation.
 *
 * Callers own both the accumulated sum and the compensation storage. Both
 * should be initialized before the first call, usually to zero.
 *
 * \param[in,out] sum Running sum updated by \p value.
 * \param[in,out] compensation Kahan compensation term updated to track the
 *                             rounding error from the addition.
 * \param[in] value Value to add to \p sum.
 */
template <typename T, typename U>
AMREX_GPU_HOST_DEVICE
void compensatedAdd (T& sum, T& compensation, U value) noexcept
{
#if defined(__clang__)
#pragma clang fp reassociate(off) contract(off)
#endif
    T const y = static_cast<T>(value) - compensation;
    T const t = sum + y;
    compensation = (t - sum) - y;
    sum = t;
}

#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
#pragma float_control(pop)
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

#if defined(__CUDACC__)
#pragma nv_diag_default 1675
#endif

}

#endif
