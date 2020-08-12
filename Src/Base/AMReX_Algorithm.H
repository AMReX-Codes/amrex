#ifndef AMREX_ALGORITHM_H_
#define AMREX_ALGORITHM_H_

#include <algorithm>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <AMReX_Dim3.H>
#include <AMReX_BLassert.H>

namespace amrex
{
    template <class T>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& min (const T& a, const T& b) noexcept
    {
        return (b < a) ? b : a;
    }

    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& min (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return min(min(a,b),c...);
    }

    template <class T>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& max (const T& a, const T& b) noexcept
    {
        return (a < b) ? b : a;
    }

    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& max (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return max(max(a,b),c...);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    Dim3 elemwiseMin (Dim3 const& a, Dim3 const& b) noexcept {
        return Dim3{amrex::min(a.x,b.x),amrex::min(a.y,b.y),amrex::min(a.z,b.z)};
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    Dim3 elemwiseMax (Dim3 const& a, Dim3 const& b) noexcept {
        return Dim3{amrex::max(a.x,b.x),amrex::max(a.y,b.y),amrex::max(a.z,b.z)};
    }

    template<typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Swap (T& t1, T& t2) noexcept
    {
        T temp = std::move(t1);
        t1 = std::move(t2);
        t2 = std::move(temp);
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr const T& Clamp (const T& v, const T& lo, const T& hi)
    {
        return (v < lo) ? lo : (hi < v) ? hi : v;
    }

    template <class T, class F>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T bisect (T lo, T hi, F f, T tol=1e-12, int max_iter=100)
    {
        AMREX_ASSERT_WITH_MESSAGE(hi > lo,
            "Error - calling bisect but lo and hi don't describe a reasonable interval.");

        T flo = f(lo);
        T fhi = f(hi);

        if (flo == T(0)) return flo;
        if (fhi == T(0)) return fhi;

        AMREX_ASSERT_WITH_MESSAGE(flo * fhi <= T(0),
            "Error - calling bisect but lo and hi don't bracket a root.");

        T mi = (lo + hi) / T(2);
        T fmi = 0.0;
        int n = 1;
        while (n <= max_iter)
        {
            if (hi - lo < tol) break;
            mi = (lo + hi) / T(2);
            fmi = f(mi);
            if (fmi == T(0)) break;
            fmi*flo < T(0) ? hi = mi : lo = mi;
            flo = f(lo);
            fhi = f(hi);
            ++n;
        }

        AMREX_ASSERT_WITH_MESSAGE(n < max_iter,
            "Error - maximum number of iterations reached in bisect.");

        return mi;
    }
}

#endif
