#ifndef AMREX_DIM3_H_
#define AMREX_DIM3_H_
#include <AMReX_Config.H>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

#include <iostream>

namespace amrex {

struct Dim3 { int x; int y; int z; };
struct XDim3 { Real x; Real y; Real z; };

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real dot_product (XDim3 const& a, XDim3 const& b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
XDim3 cross_product (XDim3 const& a, XDim3 const& b)
{
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
XDim3 operator+ (XDim3 const& a, XDim3 const& b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
XDim3 operator- (XDim3 const& a, XDim3 const& b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

template <typename T,
          std::enable_if_t <std::is_same_v<T,  Dim3> ||
                                   std::is_same_v<T, XDim3>>* = nullptr>
std::ostream& operator<< (std::ostream& os, const T& d)
{
    os << '(' << d.x << ',' << d.y << ',' << d.z << ')';
    return os;
}

}

#endif
