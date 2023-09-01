#ifndef AMREX_FUNCTIONAL_H_
#define AMREX_FUNCTIONAL_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>

// Note that since C++14, std::functional class's operator() is constexpr.
// So we don't need to use the classes here except for Minimum and Maximum.

namespace amrex {

template <typename T>
struct Plus
{
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return lhs + rhs;
    }
};

template <typename T>
struct Minus
{
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return lhs - rhs;
    }
};

template <typename T>
struct Minimum
{
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return (lhs < rhs) ? lhs : rhs;
    }
};

template <typename T>
struct Maximum
{
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return (lhs > rhs) ? lhs : rhs;
    }
};

template <typename T>
struct LogicalAnd
{
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return lhs && rhs;
    }
};

template <typename T>
struct LogicalOr
{
    constexpr T operator() (const T & lhs, const T & rhs) const
    {
        return lhs || rhs;
    }
};

}

#endif
