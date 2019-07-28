#ifndef AMREX_FUNCTIONAL_H_
#define AMREX_FUNCTIONAL_H_

#include <AMReX_GpuQualifiers.H>

namespace amrex {

template <typename T>
struct Plus
{
    constexpr T operator()(const T &lhs, const T &rhs) const
    {
        return lhs + rhs;
    }
};

}

#endif
