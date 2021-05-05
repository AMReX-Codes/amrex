#ifndef AMREX_DIM3_H_
#define AMREX_DIM3_H_
#include <AMReX_Config.H>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

#include <iostream>

namespace amrex {

struct Dim3 { int x; int y; int z; };
struct XDim3 { Real x; Real y; Real z; };

template <typename T,
          typename std::enable_if <std::is_same<T,  Dim3>::value ||
                                   std::is_same<T, XDim3>::value>::type* = nullptr>
std::ostream& operator<< (std::ostream& os, const T& d)
{
    os << '(' << d.x << ',' << d.y << ',' << d.z << ')';
    return os;
}

}

#endif
