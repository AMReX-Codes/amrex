#ifndef AMREX_ALGORITHM_H_
#define AMREX_ALGORITHM_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_Algorithm.H
 * \brief General-purpose algorithm utilities available on both host and device.
 */

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <AMReX_Dim3.H>
#include <AMReX_BLassert.H>
#include <AMReX_Math.H>

#include <algorithm>
#include <concepts>
#include <limits>
#include <type_traits>
#include <cstdint>
#include <climits>

namespace amrex
{
    /// \ingroup amrex_utilities
    //! Return the smaller value. This function was added to AMReX to
    //! support GPU before `std::min` was constexpr in C++14. `std::min` can
    //! now be used directly instead.
    template <class T>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& min (const T& a, const T& b) noexcept
    {
        return std::min(a,b);
    }

    /// \ingroup amrex_utilities
    //! Return the smallest value. This function was added to AMReX to
    //! support GPU before `std::min` was constexpr in C++14. `std::min` can
    //! now be used directly instead.
    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& min (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return min(min(a,b),c...);
    }

    /// \ingroup amrex_utilities
    //! Return the greater value. This function was added to AMReX to
    //! support GPU before `std::max` was constexpr in C++14. `std::max` can
    //! now be used directly instead.
    template <class T>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& max (const T& a, const T& b) noexcept
    {
        return std::max(a,b);
    }

    /// \ingroup amrex_utilities
    //! Return the greatest value. This function was added to AMReX to
    //! support GPU before `std::max` was constexpr in C++14. `std::max` can
    //! now be used directly instead.
    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE constexpr const T& max (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return max(max(a,b),c...);
    }

    /// \ingroup amrex_utilities
    //! Return the element-wise minimum of \p a and \p b for types like XDim3.
    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T elemwiseMin (T const& a, T const& b) noexcept {
        return T{amrex::min(a.x,b.x),amrex::min(a.y,b.y),amrex::min(a.z,b.z)};
    }

    /// \ingroup amrex_utilities
    //! Return the element-wise minimum of \p a, \p b, and the rest of the pack \p c.
    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T elemwiseMin (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return elemwiseMin(elemwiseMin(a,b),c...);
    }

    /// \ingroup amrex_utilities
    //! Return the element-wise maximum of \p a and \p b for types like XDim3.
    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T elemwiseMax (T const& a, T const& b) noexcept {
        return T{amrex::max(a.x,b.x),amrex::max(a.y,b.y),amrex::max(a.z,b.z)};
    }

    /// \ingroup amrex_utilities
    //! Return the element-wise maximum of \p a, \p b, and the remaining arguments in \p c.
    template <class T, class ... Ts>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T elemwiseMax (const T& a, const T& b, const Ts& ... c) noexcept
    {
        return elemwiseMax(elemwiseMax(a,b),c...);
    }

    /// \ingroup amrex_utilities
    //! Swap the contents of \p t1 and \p t2. `std::swap` can be used directly instead in
    //! GPU codes since C++20.
    template<typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Swap (T& t1, T& t2) noexcept
    {
        T temp = std::move(t1);
        t1 = std::move(t2);
        t2 = std::move(temp);
    }

    /// \ingroup amrex_utilities
    //! Return the reference to \p lo if \p v < \p lo; return the reference to
    //! \p hi if \p hi < \p v; otherwise return the reference to \p v. This
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

    /// \ingroup amrex_utilities
    //! Return true when the difference between \p x and \p y is small, as measured by
    //! the given \p ulp (units in the last place).
    template <std::floating_point T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool
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
     * \ingroup amrex_utilities
     * \brief Find a root of a scalar function on a bracketing interval
     * using bisection.
     *
     * Invalid bracketing and iteration failures are checked with AMReX
     * assertions when assertions are enabled.
     *
     * \tparam T floating-point type
     * \tparam F callable type of the scalar function
     * \param lo lower bound
     * \param hi upper bound
     * \param f scalar function
     * \param tol absolute tolerance. Iteration stops when `hi-lo < tol` or when
     *            `almostEqual(lo,hi)` returns true.
     * \param max_iter maximum number of bisection iterations allowed.
     * \return an approximate root found using bisection.
     */
    template <std::floating_point T, class F>
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
     * \ingroup amrex_utilities
     * \brief Find the index of the interval containing a value in a sorted
     * array.
     *
     * Find index `i` in the range [lo,hi) that d[i] <= v < d[i+1]. It is
     * assumed that the input data is sorted.
     * The elements d[lo] and d[hi] are never accessed,
     * so that d only needs to be valid in the exclusive range (lo,hi).
     * It is not checked whether the value v is out-of-bounds.
     * In that case, if v < d[lo], lo is returned, and if v >= d[hi], hi-1 is returned.
     * Note that this is different from std::lower_bound that searches
     * for the first element which is not less than `v`.
     *
     * \tparam T value type.
     * \tparam I integral index type.
     * \param d pointer to a sorted array of values.
     * \param lo inclusive lower bound of the search range.
     * \param hi exclusive upper bound of the search range.
     * \param v value to be located.
     * \return an index `i` such that `d[i] <= v < d[i+1]`, with `lo <= i < hi`.
     */
    template <typename T, std::integral I>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    I bisect (T const* d, I lo, I hi, T const& v)
    {
        while (hi - lo > 1) {
            I mid = lo + (hi - lo) / 2;
            if (v < d[mid]) {
                hi = mid;
            } else {
                lo = mid;
            }
        };
        return lo;
    }

    /**
     * \ingroup amrex_utilities
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
     * \ingroup amrex_utilities
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
     * \ingroup amrex_utilities
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
    template<typename ItType, std::floating_point ValType>
    requires (std::floating_point<typename std::iterator_traits<ItType>::value_type>)
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
     * \ingroup amrex_utilities
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
    template<typename ItType, std::floating_point ValType>
    requires (std::floating_point<typename std::iterator_traits<ItType>::value_type>)
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
}

#endif
