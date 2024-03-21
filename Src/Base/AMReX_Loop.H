#ifndef AMREX_LOOP_H_
#define AMREX_LOOP_H_
#include <AMReX_Config.H>

#include <AMReX_Box.H>

namespace amrex {

template <class F>
AMREX_GPU_HOST_DEVICE
void Loop (Dim3 lo, Dim3 hi, F const& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_GPU_HOST_DEVICE
void Loop (Dim3 lo, Dim3 hi, int ncomp, F const& f) noexcept
{
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
AMREX_GPU_HOST_DEVICE
void LoopConcurrent (Dim3 lo, Dim3 hi, F const& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_GPU_HOST_DEVICE
void LoopConcurrent (Dim3 lo, Dim3 hi, int ncomp, F const& f) noexcept
{
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
AMREX_GPU_HOST_DEVICE
void Loop (Box const& bx, F const& f) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_GPU_HOST_DEVICE
void Loop (Box const& bx, int ncomp, F const& f) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
AMREX_GPU_HOST_DEVICE
void LoopConcurrent (Box const& bx, F const& f) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_GPU_HOST_DEVICE
void LoopConcurrent (Box const& bx, int ncomp, F const& f) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

// The functions above are __host__ __device__ functions.  If f is not a
// __host__ __device__ function, we will get warning about calling __host__
// function from a __host__ __device__ function.  This is ugly.  To get rid
// of the warning, we have to use the functions below for those situations.

template <class F>
void LoopOnCpu (Dim3 lo, Dim3 hi, F const& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
void LoopOnCpu (Dim3 lo, Dim3 hi, int ncomp, F const& f) noexcept
{
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
void LoopConcurrentOnCpu (Dim3 lo, Dim3 hi, F const& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
void LoopConcurrentOnCpu (Dim3 lo, Dim3 hi, int ncomp, F const& f) noexcept
{
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
void LoopOnCpu (Box const& bx, F const& f) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
void LoopOnCpu (Box const& bx, int ncomp, F const& f) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
void LoopConcurrentOnCpu (Box const& bx, F const& f) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
void LoopConcurrentOnCpu (Box const& bx, int ncomp, F const& f) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

// Implementation of "constexpr for" based on
// https://artificial-mind.net/blog/2020/10/31/constexpr-for
//
// Approximates what one would get from a compile-time
// unrolling of the loop
// for (int i = 0; i < N; ++i) {
//    f(i);
// }
//
// The mechanism is recursive: we evaluate f(i) at the current
// i and then call the for loop at i+1. f() is a lambda function
// that provides the body of the loop and takes only an integer
// i as its argument.

template<auto I, auto N, class F>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr void constexpr_for (F const& f)
{
    if constexpr (I < N) {
        f(std::integral_constant<decltype(I), I>());
        constexpr_for<I+1, N>(f);
    }
}

#include <AMReX_Loop.nolint.H>

}

#endif
