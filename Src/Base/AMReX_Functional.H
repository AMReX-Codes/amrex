#ifndef AMREX_FUNCTIONAL_H_
#define AMREX_FUNCTIONAL_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_Functional.H
 * \brief Function object (functor) types for use with reductions and transforms.
 */

#include <AMReX_GpuQualifiers.H>

// Note that since C++14, std::functional class's operator() is constexpr.
// So we don't need to use the classes here except for Minimum and Maximum.

namespace amrex {

/**
 * \ingroup amrex_utilities
 * \brief Function object that returns the sum of two values.
 */
template <typename T>
struct Plus
{
    //! Return lhs + rhs.
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return lhs + rhs;
    }
};

/**
 * \ingroup amrex_utilities
 * \brief Function object that returns the difference of two values.
 */
template <typename T>
struct Minus
{
    //! Return lhs - rhs.
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return lhs - rhs;
    }
};

/**
 * \ingroup amrex_utilities
 * \brief Function object that returns the smaller of two values.
 */
template <typename T>
struct Minimum
{
    //! Return the smaller of \p lhs and \p rhs.
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return (lhs < rhs) ? lhs : rhs;
    }
};

/**
 * \ingroup amrex_utilities
 * \brief Function object that returns the larger of two values.
 */
template <typename T>
struct Maximum
{
    //! Return the larger of \p lhs and \p rhs.
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return (lhs > rhs) ? lhs : rhs;
    }
};

/**
 * \ingroup amrex_utilities
 * \brief Function object that returns the logical AND of two values.
 */
template <typename T>
struct LogicalAnd
{
    //! Return lhs && rhs.
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return lhs && rhs;
    }
};

/**
 * \ingroup amrex_utilities
 * \brief Function object that returns the logical OR of two values.
 */
template <typename T>
struct LogicalOr
{
    //! Return lhs || rhs.
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return lhs || rhs;
    }
};

/**
 * \ingroup amrex_utilities
 * \brief Function object that returns the product of two values.
 */
template <typename T>
struct Multiplies
{
    //! Return lhs * rhs.
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return lhs * rhs;
    }
};

/**
 * \ingroup amrex_utilities
 * \brief Function object that returns the quotient of two values.
 */
template <typename T>
struct Divides
{
    //! Return lhs / rhs.
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return lhs / rhs;
    }
};

}

#endif
