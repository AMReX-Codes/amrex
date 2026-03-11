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
    //! Return the smaller value. This function was added to AMReX to
    //! support GPU before `std::min` was constexpr in C++14. `std::min` can
    //! now be used directly instead.
    template <class T>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& min (const T& a, const T& b) noexcept
    {
        return std::min(a,b);
    }

    //! Return the smallest value. This function was added to AMReX to
    //! support GPU before `std::min` was constexpr in C++14. `std::min` can
    //! now be used directly instead.
    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& min (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return min(min(a,b),c...);
    }

    //! Return the greater value. This function was added to AMReX to
    //! support GPU before `std::max` was constexpr in C++14. `std::max` can
    //! now be used directly instead.
    template <class T>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& max (const T& a, const T& b) noexcept
    {
        return std::max(a,b);
    }

    //! Return the greatest value. This function was added to AMReX to
    //! support GPU before `std::max` was constexpr in C++14. `std::max` can
    //! now be used directly instead.
    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& max (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return max(max(a,b),c...);
    }

    //! Return the element-wise minimum of the given values for types like XDim3.
    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T elemwiseMin (T const& a, T const& b) noexcept {
        return T{amrex::min(a.x,b.x),amrex::min(a.y,b.y),amrex::min(a.z,b.z)};
    }

    //! Return the element-wise minimum of the given values for types like XDim3.
    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T elemwiseMin (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return elemwiseMin(elemwiseMin(a,b),c...);
    }

    //! Return the element-wise maximum of the given values for types like XDim3.
    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T elemwiseMax (T const& a, T const& b) noexcept {
        return T{amrex::max(a.x,b.x),amrex::max(a.y,b.y),amrex::max(a.z,b.z)};
    }

    //! Return the element-wise maximum of the given values for types like XDim3.
    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T elemwiseMax (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return elemwiseMax(elemwiseMax(a,b),c...);
    }

    //! Swap the given values. `std::swap` can be used directly instead in
    //! GPU codes since C++20.
    template<typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Swap (T& t1, T& t2) noexcept
    {
        T temp = std::move(t1);
        t1 = std::move(t2);
        t2 = std::move(temp);
    }

    //! Return the reference to `lo` if `v < lo`; return the reference to
    //! `hi` if `hi < v`; otherwise return the reference to `v. This
    //! function was added to AMReX before switching to C++17. `std::clamp`
    //! can now be used directly instead.
    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr const T& Clamp (const T& v, const T& lo, const T& hi)
    {
        if (v < lo) {
            return lo; // NOLINT(bugprone-return-const-ref-from-parameter)
        } else if (hi < v) {
            return hi; // NOLINT(bugprone-return-const-ref-from-parameter)
        } else {
            return v; // NOLINT(bugprone-return-const-ref-from-parameter)
        }
    }

    //! Return if the difference between two values are small as measured by
    //! the given `ulp` (units in the last place).
    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<std::is_floating_point_v<T>,bool>
    almostEqual (T x, T y, int ulp = 2)
    {
        // Reference: https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
        // The machine epsilon has to be scaled to the magnitude of the values used
        // and multiplied by the desired precision in ULPs (units in the last place)
        return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
        // unless the result is subnormal
        || std::abs(x-y) < std::numeric_limits<T>::min();
    }

    /**
     * \brief Find a root of a scalar function on a bracketing interval
     * using bisection.
     *
     * It's a runtime error if the root finding fails.
     *
     * \tparam T floating-point type
     * \tparam F callable type of the scalar function
     * \param lo lower bound
     * \param hi upper bound
     * \param f scalar function
     * \param tol absolute tolerance. Iteration stops when `hi-lo < tol` or when
     *            `almostEqual(lo,hi) returns true.
     * \param max_iter maximum number of bisection iterations allowed.
     * \return an approximate root found using bisection.
     */
    template <class T, class F,
              std::enable_if_t<std::is_floating_point_v<T>,int>FOO = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T bisect (T lo, T hi, F f, T tol=1e-12, int max_iter=100)
    {
        AMREX_ASSERT_WITH_MESSAGE(hi > lo,
            "Error - calling bisect but lo and hi don't describe a reasonable interval.");

        T flo = f(lo);
        T fhi = f(hi);

        if (flo == T(0)) { return lo; }
        if (fhi == T(0)) { return hi; }

        AMREX_ASSERT_WITH_MESSAGE(flo * fhi <= T(0),
            "Error - calling bisect but lo and hi don't bracket a root.");

        T mi = (lo + hi) / T(2);
        T fmi = T(0);
        int n = 1;
        while (n <= max_iter)
        {
            if (hi - lo < tol || almostEqual(lo,hi)) { break; }
            mi = (lo + hi) / T(2);
            fmi = f(mi);
            if (fmi == T(0)) { break; }
            fmi*flo < T(0) ? hi = mi : lo = mi;
            flo = f(lo);
            fhi = f(hi);
            ++n;
        }

        AMREX_ASSERT_WITH_MESSAGE(n <= max_iter,
            "Error - maximum number of iterations reached in bisect.");

        return mi;
    }

    /**
     * \brief Find the index of the interval containing a value in a sorted
     * array.
     *
     * Find index `i` in the range [lo,hi) that d[i] <= v < d[i+1].  It is
     * assumed that the input data are sorted and d[lo] <= v < d[hi].  Note
     * that this is different from std::lower_bound that searches for the
     * first element which is not less than `v`.
     *
     * \tparam T value type.
     * \tparam I integral index type.
     * \param d pointer to a sorted array of values.
     * \param lo inclusive lower bound of the search range.
     * \param hi exclusive upper bound of the search range.
     * \param v value to be located.
     * \return an index `i` such that `d[i] <= v < d[i+1]`.
     */
    template <typename T, typename I,
              std::enable_if_t<std::is_integral_v<I>,int> = 0>
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

    /**
     * \brief Return an iterator to the first element greater than a given
     * value.
     *
     * This function is an implementation of `std::upper_bound` that works
     * on both host and device.
     *
     * \tparam ItType iterator type.
     * \tparam ValType value type.
     * \param first inclusive lower bound of the search range.
     * \param last exclusive upper bound of the search range.
     * \param val value to compare the elements to.
     * \return an iterator pointing to the first element greater than `val`.
     */
    template<typename ItType, typename ValType>
    AMREX_GPU_HOST_DEVICE
    ItType upper_bound (ItType first, ItType last, const ValType& val)
    {
        AMREX_IF_ON_DEVICE((
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
        ))
        AMREX_IF_ON_HOST((
            return std::upper_bound(first, last, val);
        ))
    }

    /**
     * \brief Return an iterator to the first element not less than a given
     * value.
     *
     * This function is an implementation of `std::lower_bound` that works
     * on both host and device.
     *
     * \tparam ItType iterator type.
     * \tparam ValType value type.
     * \param first inclusive lower bound of the search range.
     * \param last exclusive upper bound of the search range.
     * \param val value to compare the elements to.
     * \return an iterator pointing to the first element not less than `val`.
     */
    template<typename ItType, typename ValType>
    AMREX_GPU_HOST_DEVICE
    ItType lower_bound (ItType first, ItType last, const ValType& val)
    {
        AMREX_IF_ON_DEVICE((
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
        ))
        AMREX_IF_ON_HOST((
            return std::lower_bound(first, last, val);
        ))
    }

    /**
     * \brief Fill a range with linearly spaced values over a closed interval.
     *
     * This function assigns linearly spaced floating-point values to the
     * range `[first, last)`, starting at `start` and ending at `stop`. The
     * value range is inclusive at both ends such that the first element is
     * set to `start` and the last element to `stop` exactly. Note that this
     * function does nothing when the range contains fewer than two elements
     * (i.e., `last-first < 2`).
     *
     * \tparam ItType iterator type.
     * \tparam ValType floating-point value type.
     * \param first pointing to the first element of the output range.
     * \param last pointing one past the last element of the output range.
     * \param start start value.
     * \param stop stop value.
     */
    template<typename ItType, typename ValType,
        std::enable_if_t<
            std::is_floating_point_v<typename std::iterator_traits<ItType>::value_type> &&
            std::is_floating_point_v<ValType>,
        int> = 0>
    AMREX_GPU_HOST_DEVICE
    void linspace (ItType first, const ItType& last, const ValType& start, const ValType& stop)
    {
        const std::ptrdiff_t count = last-first;
        if (count >= 2){
            const auto delta = (stop - start)/(count - 1);
            for (std::ptrdiff_t i = 0; i < count-1; ++i){
                *(first++) = start + i*delta;
            }
            *first = stop;
        }
    }

    /**
     * \brief Fill a range with logarithmically spaced values over a closed interval.
     *
     * This function assigns logarithmically spaced floating-point values to
     * the range `[first, last)`, starting at `base^start` and ending at
     * `base^stop`. The value range is inclusive at both ends such that the
     * first element is set to `base^start` and the last element to
     * `base^stop` exactly. Note that this function does nothing when the
     * range contains fewer than two elements (i.e., `last-first < 2`).
     *
     * \tparam ItType iterator type.
     * \tparam ValType floating-point value type.
     * \param first pointing to the first element of the output range.
     * \param last pointing one past the last element of the output range.
     * \param start start value.
     * \param stop stop value.
     * \param base base of the exponential.
     */
    template<typename ItType, typename ValType,
        std::enable_if_t<
            std::is_floating_point_v<typename std::iterator_traits<ItType>::value_type> &&
            std::is_floating_point_v<ValType>,
        int> = 0>
    AMREX_GPU_HOST_DEVICE
    void logspace (ItType first, const ItType& last,
        const ValType& start, const ValType& stop, const ValType& base)
    {
        const std::ptrdiff_t count = last-first;
        if (count >= 2){
            const auto delta = (stop - start)/(count - 1);
            for (std::ptrdiff_t i = 0; i < count-1; ++i){
                *(first++) = std::pow(base, start + i*delta);
            }
            *first = std::pow(base, stop);
        }
    }

/// \cond DOXYGEN_IGNORE
namespace detail {

struct clzll_tag {};
struct clzl_tag : clzll_tag {};
struct clz_tag : clzl_tag {};

// in gcc and clang, there are three versions of __builtin_clz taking unsigned int,
// unsigned long, and unsigned long long inputs. Because the sizes of these data types
// vary on different platforms, we work with fixed-width integer types.
// these tags and overloads select the smallest version of __builtin_clz that will hold the input type
template <typename T, typename = std::enable_if_t<sizeof(T) <= sizeof(unsigned int)>>
AMREX_FORCE_INLINE
int builtin_clz_wrapper (clz_tag, T x) noexcept
{
    return static_cast<int>(__builtin_clz(x) - (sizeof(unsigned int) * CHAR_BIT - sizeof(T) * CHAR_BIT));
}

template <typename T, typename = std::enable_if_t<sizeof(T) <= sizeof(unsigned long)>>
AMREX_FORCE_INLINE
int builtin_clz_wrapper (clzl_tag, T x) noexcept
{
    return static_cast<int>(__builtin_clzl(x) - (sizeof(unsigned long) * CHAR_BIT - sizeof(T) * CHAR_BIT));
}

template <typename T, typename = std::enable_if_t<sizeof(T) <= sizeof(unsigned long long)>>
AMREX_FORCE_INLINE
int builtin_clz_wrapper (clzll_tag, T x) noexcept
{
    return static_cast<int>(__builtin_clzll(x) - (sizeof(unsigned long long) * CHAR_BIT - sizeof(T) * CHAR_BIT));
}

}
/// \endcond

//! Return the number of leading zeros of the given integer.
template <class T, std::enable_if_t<std::is_same_v<std::decay_t<T>,std::uint8_t>  ||
                                    std::is_same_v<std::decay_t<T>,std::uint16_t> ||
                                    std::is_same_v<std::decay_t<T>,std::uint32_t> ||
                                    std::is_same_v<std::decay_t<T>,std::uint64_t>, int> = 0>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clz (T x) noexcept;

/// \cond DOXYGEN_IGNORE

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clz_generic (std::uint8_t x) noexcept
{
#if !defined(__NVCOMPILER)
    static constexpr int clz_lookup[16] = { 4, 3, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
#else
    constexpr int clz_lookup[16] = { 4, 3, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
#endif
    auto upper = x >> 4;
    auto lower = x & 0xF;
    return upper ? clz_lookup[upper] : 4 + clz_lookup[lower];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clz_generic (std::uint16_t x) noexcept
{
    auto upper = std::uint8_t(x >> 8);
    auto lower = std::uint8_t(x & 0xFF);
    return upper ? clz(upper) : 8 + clz(lower);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clz_generic (std::uint32_t x) noexcept
{
    auto upper = std::uint16_t(x >> 16);
    auto lower = std::uint16_t(x & 0xFFFF);
    return upper ? clz(upper) : 16 + clz(lower);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clz_generic (std::uint64_t x) noexcept
{
    auto upper = std::uint32_t(x >> 32);
    auto lower = std::uint32_t(x & 0xFFFFFFFF);
    return upper ? clz(upper) : 32 + clz(lower);
}

/// \endcond

#if defined AMREX_USE_CUDA

/// \cond DOXYGEN_IGNORE
namespace detail {
    // likewise with CUDA, there are __clz functions that take (signed) int and long long int
    template <typename T, typename = std::enable_if_t<sizeof(T) <= sizeof(int)> >
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    int clz_wrapper (clz_tag, T x) noexcept
    {
        return __clz((int) x) - (sizeof(int) * CHAR_BIT - sizeof(T) * CHAR_BIT);
    }

    template <typename T, typename = std::enable_if_t<sizeof(T) <= sizeof(long long int)> >
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    int clz_wrapper (clzll_tag, T x) noexcept
    {
        return __clzll((long long int) x) - (sizeof(long long int) * CHAR_BIT - sizeof(T) * CHAR_BIT);
    }
}
/// \endcond

template <class T, std::enable_if_t<std::is_same_v<std::decay_t<T>,std::uint8_t>  ||
                                    std::is_same_v<std::decay_t<T>,std::uint16_t> ||
                                    std::is_same_v<std::decay_t<T>,std::uint32_t> ||
                                    std::is_same_v<std::decay_t<T>,std::uint64_t>, int> >
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clz (T x) noexcept
{
    AMREX_IF_ON_DEVICE((return detail::clz_wrapper(detail::clz_tag{}, x);))
#if AMREX_HAS_BUILTIN_CLZ
    AMREX_IF_ON_HOST((return detail::builtin_clz_wrapper(detail::clz_tag{}, x);))
#else
    AMREX_IF_ON_HOST((return clz_generic(x);))
#endif
}

#else // !defined AMREX_USE_CUDA

template <class T, std::enable_if_t<std::is_same_v<std::decay_t<T>,std::uint8_t>  ||
                                    std::is_same_v<std::decay_t<T>,std::uint16_t> ||
                                    std::is_same_v<std::decay_t<T>,std::uint32_t> ||
                                    std::is_same_v<std::decay_t<T>,std::uint64_t>, int> >
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clz (T x) noexcept
{
#if (!AMREX_DEVICE_COMPILE && AMREX_HAS_BUILTIN_CLZ)
    return detail::builtin_clz_wrapper(detail::clz_tag{}, x);
#else
    return clz_generic(x);
#endif
}

#endif // defined AMREX_USE_CUDA

}

#endif
