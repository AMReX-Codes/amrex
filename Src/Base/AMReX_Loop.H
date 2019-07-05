#ifndef AMREX_LOOP_H_
#define AMREX_LOOP_H_

#include <AMReX_Box.H>

namespace amrex {

template <class F>
AMREX_GPU_HOST_DEVICE
void Loop (Dim3 lo, Dim3 hi, F&& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_GPU_HOST_DEVICE
void Loop (Dim3 lo, Dim3 hi, int ncomp, F&& f) noexcept
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
void LoopConcurrent (Dim3 lo, Dim3 hi, F&& f) noexcept
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
void LoopConcurrent (Dim3 lo, Dim3 hi, int ncomp, F&& f) noexcept
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
void Loop (Box const& bx, F&& f) noexcept
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
void Loop (Box const& bx, int ncomp, F&& f) noexcept
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
void LoopConcurrent (Box const& bx, F&& f) noexcept
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
void LoopConcurrent (Box const& bx, int ncomp, F&& f) noexcept
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
void LoopOnCpu (Dim3 lo, Dim3 hi, F&& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
void LoopOnCpu (Dim3 lo, Dim3 hi, int ncomp, F&& f) noexcept
{
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
void LoopConcurrentOnCpu (Dim3 lo, Dim3 hi, F&& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
void LoopConcurrentOnCpu (Dim3 lo, Dim3 hi, int ncomp, F&& f) noexcept
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
void LoopOnCpu (Box const& bx, F&& f) noexcept
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
void LoopOnCpu (Box const& bx, int ncomp, F&& f) noexcept
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
void LoopConcurrentOnCpu (Box const& bx, F&& f) noexcept
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
void LoopConcurrentOnCpu (Box const& bx, int ncomp, F&& f) noexcept
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

}

#endif
