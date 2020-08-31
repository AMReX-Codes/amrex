#ifndef AMREX_DIM3_H_
#define AMREX_DIM3_H_

#include <iostream>
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace amrex {

struct Dim3 { int x; int y; int z; };
struct XDim3 { Real x; Real y; Real z; };

inline std::ostream& operator<< (std::ostream& os, const Dim3& d) {
    os << '(' << d.x << ',' << d.y << ',' << d.z << ')';
    return os;
}

}

#endif
