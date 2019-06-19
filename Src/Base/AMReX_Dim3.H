#ifndef AMREX_DIM3_H_
#define AMREX_DIM3_H_

#include <iostream>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>

namespace amrex {

struct Dim3 { int x; int y; int z; };
struct XDim3 { Real x; Real y; Real z; };

inline std::ostream& operator<< (std::ostream& os, const Dim3& d) {
    os << '(' << d.x << ',' << d.y << ',' << d.z << ')';
    return os;
}

inline XDim3 makeXDim3 (const Array<Real,AMREX_SPACEDIM>& a) noexcept
{
#if (AMREX_SPACEDIM == 1)
    return XDim3{a[0], 0., 0.};
#elif (AMREX_SPACEDIM == 2)
    return XDim3{a[0], a[1], 0.};
#else
    return XDim3{a[0], a[1], a[2]};
#endif
}

}

#endif
