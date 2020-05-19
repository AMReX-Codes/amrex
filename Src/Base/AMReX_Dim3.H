#ifndef AMREX_DIM3_H_
#define AMREX_DIM3_H_

#include <iostream>
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace amrex {

#ifdef AMREX_USE_HIP
struct Dim3 { int x; int y; int z;
    AMREX_GPU_HOST_DEVICE
    constexpr Dim3 () : x(0), y(0), z(0) {}
    AMREX_GPU_HOST_DEVICE
    constexpr Dim3 (int x_, int y_, int z_) // xxxxx HIP: todo
        : x(x_), y(y_), z(z_) {}
    AMREX_GPU_HOST_DEVICE
    constexpr Dim3 (Dim3 const& rhs)
        : x(rhs.x), y(rhs.y), z(rhs.z) {}
};
struct XDim3 { Real x; Real y; Real z;
    AMREX_GPU_HOST_DEVICE
    constexpr XDim3 () : x(0._rt), y(0._rt), z(0._rt) {}
    AMREX_GPU_HOST_DEVICE
    constexpr XDim3 (Real x_, Real y_, Real z_) // xxxxx HIP: todo
        : x(x_), y(y_), z(z_) {}
    AMREX_GPU_HOST_DEVICE
    constexpr XDim3 (XDim3 const& rhs)
        : x(rhs.x), y(rhs.y), z(rhs.z) {}
};
#else
struct Dim3 { int x; int y; int z; };
struct XDim3 { Real x; Real y; Real z; };
#endif

inline std::ostream& operator<< (std::ostream& os, const Dim3& d) {
    os << '(' << d.x << ',' << d.y << ',' << d.z << ')';
    return os;
}

}

#endif
