#ifndef AMREX_COORDSYS_1D_C_H_
#define AMREX_COORDSYS_1D_C_H_

#include <AMReX_Gpu.H>
#include <AMReX_Array.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

namespace amrex {

AMREX_GPU_HOST_DEVICE
inline
void
amrex_setvol (Box const& bx, FArrayBox& vol,
                  GpuArray<Real,1> const& offset,
                  GpuArray<Real,1> const& dx, const int coord)
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = vol.view(lo);

    if (coord == 0) // Cartesian
    {
        Real dv = dx[0];
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            dp(i,0,0) = dv;
        }
    }
    else if (coord == 1)
    {
        const Real pi = 3.1415926535897932;
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            Real ri = offset[0] + dx[0]*(i+lo.x);
            Real ro = ri + dx[0];
            Real v = pi*(ro-ri)*(ro + ri);
            dp(i,0,0) = std::abs(v);
        }
    }
    else
    {
        const Real pi = 3.1415926535897932;
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            Real ri = offset[0] + dx[0]*(i+lo.x);
            Real ro = ri + dx[0];
            Real v = ((4./3.)*pi)*(ro-ri)*(ro*ro+ro*ri+ri*ri);
            dp(i,0,0) = std::abs(v);
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void
amrex_setarea (Box const& bx, FArrayBox& area,
                   GpuArray<Real,1> const& offset,
                   GpuArray<Real,1> const& dx, const int /*dir*/, const int coord)
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = area.view(lo);

    if (coord == 0)
    {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            dp(i,0,0) = 1.0;
        }
    }
    else if (coord == 1)
    {
        const Real pi = 3.1415926535897932;
        AMREX_PRAGMA_SIMD 
        for (int i = 0; i < len.x; ++i) {
            Real ri = offset[0] + dx[0]*(i+lo.x);
            Real a = (2.*pi)*ri;
            dp(i,0,0) = std::abs(a);
        }
    }
    else
    {
        const Real pi = 3.1415926535897932;
        AMREX_PRAGMA_SIMD 
        for (int i = 0; i < len.x; ++i) {
            Real ri = offset[0] + dx[0]*(i+lo.x);
            Real a = (4.0*pi)*ri*ri;
            dp(i,0,0) = std::abs(a);
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void
amrex_setdloga (Box const& bx, FArrayBox& dloga,
                    GpuArray<Real,1> const& offset,
                    GpuArray<Real,1> const& dx, const int /*dir*/, const int coord)
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = dloga.view(lo);

    if (coord == 0)
    {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            dp(i,0,0) = 0.0;
        }
    }
    else if (coord == 1)
    {
        AMREX_PRAGMA_SIMD 
        for (int i = 0; i < len.x; ++i) {
            Real rc = offset[0] + dx[0]*(i+lo.x+0.5);
            dp(i,0,0) = 1.0/rc;
        }
    }
    else
    {
        AMREX_PRAGMA_SIMD 
        for (int i = 0; i < len.x; ++i) {
            Real ri = offset[0] + dx[0]*(i+lo.x+0.5);
            dp(i,0,0) = 2.0/ri;
        }
    }
}

}

#endif
