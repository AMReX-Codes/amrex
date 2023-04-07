#ifndef AMREX_ALGORITHM_H_
#define AMREX_ALGORITHM_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <AMReX_Dim3.H>
#include <AMReX_BLassert.H>
#include <AMReX_Math.H>

#include <algorithm>
#include <limits>
#include <type_traits>
#include <cstdint>
#include <climits>

namespace amrex
{
    template <class T>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& min (const T& a, const T& b) noexcept
    {
        return std::min(a,b);
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
        return std::max(a,b);
    }

    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& max (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return max(max(a,b),c...);
    }

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T elemwiseMin (T const& a, T const& b) noexcept {
        return T{amrex::min(a.x,b.x),amrex::min(a.y,b.y),amrex::min(a.z,b.z)};
    }

    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T elemwiseMin (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return elemwiseMin(elemwiseMin(a,b),c...);
    }

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T elemwiseMax (T const& a, T const& b) noexcept {
        return T{amrex::max(a.x,b.x),amrex::max(a.y,b.y),amrex::max(a.z,b.z)};
    }

    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T elemwiseMax (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return elemwiseMax(elemwiseMax(a,b),c...);
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

    // Reference: https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    typename std::enable_if<std::is_floating_point<T>::value,bool>::type
    almostEqual (T x, T y, int ulp = 2)
    {
        // the machine epsilon has to be scaled to the magnitude of the values used
        // and multiplied by the desired precision in ULPs (units in the last place)
        return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
        // unless the result is subnormal
        || std::abs(x-y) < std::numeric_limits<T>::min();
    }

    template <class T, class F,
              typename std::enable_if<std::is_floating_point<T>::value,int>::type FOO = 0>
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
        T fmi = T(0);
        int n = 1;
        while (n <= max_iter)
        {
            if (hi - lo < tol || almostEqual(lo,hi)) break;
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

    // Find I in the range [lo,hi) that T[I] <= v < T[I+1].
    // It is assumed that the input data are sorted and T[lo] <= v < T[hi].
    // Note that this is different from std::lower_bound.
    template <typename T, typename I,
              typename std::enable_if<std::is_integral<I>::value,int>::type = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    I bisect (T const* d, I lo, I hi, T const& v) {
        while (lo <= hi) {
            int mid = lo + (hi-lo)/2;
            if (v >= d[mid] && v < d[mid+1]) {
                return mid;
            } else if (v < d[mid]) {
                hi = mid-1;
            } else {
                lo = mid+1;
            }
        };
        return hi;
    }

    template<typename ItType, typename ValType>
    AMREX_GPU_HOST_DEVICE
    ItType upper_bound (ItType first, ItType last, const ValType& val)
    {
#if AMREX_DEVICE_COMPILE
        std::ptrdiff_t count = last-first;
        while(count>0){
            auto it = first;
            const auto step = count/2;
            it += step;
            if (!(val < *it)){
                first = ++it;
                count -= step + 1;
            }
            else{
                count = step;
            }
        }

        return first;
#else
        return std::upper_bound(first, last, val);
#endif
    }

    template<typename ItType, typename ValType>
    AMREX_GPU_HOST_DEVICE
    ItType lower_bound (ItType first, ItType last, const ValType& val)
    {
#ifdef AMREX_DEVICE_COMPILE
        std::ptrdiff_t count = last-first;
        while(count>0)
        {
            auto it = first;
            const auto step = count/2;
            it += step;
            if (*it < val){
                first = ++it;
                count -= step + 1;
            }
            else{
                count = step;
            }
        }

        return first;
#else
        return std::lower_bound(first, last, val);
#endif
    }

namespace detail {

struct clzll_tag {};
struct clzl_tag : clzll_tag {};
struct clz_tag : clzl_tag {};

// in gcc and clang, there are three versions of __builtin_clz taking unsigned int,
// unsigned long, and unsigned long long inputs. Because the sizes of these data types
// vary on different platforms, we work with fixed-width integer types.
// these tags and overloads select the smallest version of __builtin_clz that will hold the input type
template <typename T, typename = typename std::enable_if<sizeof(T) <= sizeof(unsigned int)>::type>
AMREX_FORCE_INLINE
int builtin_clz_wrapper (clz_tag, T x) noexcept
{
    return static_cast<int>(__builtin_clz(x) - (sizeof(unsigned int) * CHAR_BIT - sizeof(T) * CHAR_BIT));
}

template <typename T, typename = typename std::enable_if<sizeof(T) <= sizeof(unsigned long)>::type>
AMREX_FORCE_INLINE
int builtin_clz_wrapper (clzl_tag, T x) noexcept
{
    return static_cast<int>(__builtin_clzl(x) - (sizeof(unsigned long) * CHAR_BIT - sizeof(T) * CHAR_BIT));
}

template <typename T, typename = typename std::enable_if<sizeof(T) <= sizeof(unsigned long long)>::type>
AMREX_FORCE_INLINE
int builtin_clz_wrapper (clzll_tag, T x) noexcept
{
    return static_cast<int>(__builtin_clzll(x) - (sizeof(unsigned long long) * CHAR_BIT - sizeof(T) * CHAR_BIT));
}

#ifdef AMREX_USE_CUDA

// likewise with CUDA, there are __clz functions that take (signed) int and long long int
template <typename T, typename = typename std::enable_if<sizeof(T) <= sizeof(int)>::type>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int clz_wrapper (clz_tag, T x) noexcept
{
    return __clz((int) x) - (sizeof(int) * CHAR_BIT - sizeof(T) * CHAR_BIT);
}

template <typename T, typename = typename std::enable_if<sizeof(T) <= sizeof(long long int)>::type>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int clz_wrapper (clzll_tag, T x) noexcept
{
    return __clzll((long long int) x) - (sizeof(long long int) * CHAR_BIT - sizeof(T) * CHAR_BIT);
}
#endif

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clz (std::uint8_t x) noexcept
{
#if (AMREX_DEVICE_COMPILE && defined(AMREX_USE_CUDA)) // all supported cuda versions have __clz
    return detail::clz_wrapper(detail::clz_tag{}, x);
#elif (!AMREX_DEVICE_COMPILE && AMREX_HAS_BUILTIN_CLZ)
    return detail::builtin_clz_wrapper(detail::clz_tag{}, x);
#else
    static constexpr int clz_lookup[16] = { 4, 3, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
    auto upper = x >> 4;
    auto lower = x & 0xF;
    return upper ? clz_lookup[upper] : 4 + clz_lookup[lower];
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clz (std::uint16_t x) noexcept
{
#if (AMREX_DEVICE_COMPILE && defined(AMREX_USE_CUDA)) // all supported cuda versions have __clz
    return detail::clz_wrapper(detail::clz_tag{}, x);
#elif (!AMREX_DEVICE_COMPILE && AMREX_HAS_BUILTIN_CLZ)
    return detail::builtin_clz_wrapper(detail::clz_tag{}, x);
#else
    auto upper = std::uint8_t(x >> 8);
    auto lower = std::uint8_t(x & 0xFF);
    return upper ? clz(upper) : 8 + clz(lower);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clz (std::uint32_t x) noexcept
{
#if (AMREX_DEVICE_COMPILE && defined(AMREX_USE_CUDA)) // all supported cuda versions have __clz
    return detail::clz_wrapper(detail::clz_tag{}, x);
#elif (!AMREX_DEVICE_COMPILE && AMREX_HAS_BUILTIN_CLZ)
    return detail::builtin_clz_wrapper(detail::clz_tag{}, x);
#else
    auto upper = std::uint16_t(x >> 16);
    auto lower = std::uint16_t(x & 0xFFFF);
    return upper ? clz(upper) : 16 + clz(lower);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clz (std::uint64_t x) noexcept
{
#if (AMREX_DEVICE_COMPILE && defined(AMREX_USE_CUDA)) // all supported cuda versions have __clz
    return detail::clz_wrapper(detail::clz_tag{}, x);
#elif (!AMREX_DEVICE_COMPILE && AMREX_HAS_BUILTIN_CLZ)
    return detail::builtin_clz_wrapper(detail::clz_tag{}, x);
#else
    auto upper = std::uint32_t(x >> 32);
    auto lower = std::uint32_t(x & 0xFFFFFFFF);
    return upper ? clz(upper) : 32 + clz(lower);
#endif
}

}

#endif
