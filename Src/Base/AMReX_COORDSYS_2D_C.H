#ifndef AMREX_COORDSYS_2D_C_H_
#define AMREX_COORDSYS_2D_C_H_

#include <AMReX_Gpu.H>
#include <AMReX_Array.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

namespace amrex {

AMREX_GPU_HOST_DEVICE
inline
void
amrex_setvol (Box const& bx, Array4<Real> const& vol,
              GpuArray<Real,2> const& offset,
              GpuArray<Real,2> const& dx, const int coord) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    if (coord == 0) // Cartesian
    {
        Real dv = dx[0]*dx[1];
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                vol(i,j,0) = dv;
            }
        }
    }
    else if (coord == 1) // r-z
    {
        const Real pi = 3.1415926535897932;
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real ri = offset[0] + dx[0]*(i);
                Real ro = ri + dx[0];
                Real v = pi*dx[1]*dx[0]*(ro + ri);
                vol(i,j,0) = std::abs(v);
            }
        }
    }
    else // r-theta
    {
        const Real pi = 3.1415926535897932;
        for     (int j = lo.y; j <= hi.y; ++j) {
            Real ti = offset[1] + dx[1]*(j);
            Real to = ti + dx[1];
            Real tmp = (2.*pi)*(std::cos(ti)-std::cos(to))/3.0;
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real ri = offset[0] + dx[0]*(i);
                Real ro = ri + dx[0];
                Real v = tmp*(ro-ri)*(ro*ro+ro*ri+ri*ri);
                vol(i,j,0) = std::abs(v);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void
amrex_setarea (Box const& bx, Array4<Real> const& area,
               GpuArray<Real,2> const& offset,
               GpuArray<Real,2> const& dx, const int dir, const int coord) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    if (coord == 0)
    {
        Real a = (dir == 0) ? dx[1] : dx[0];
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                area(i,j,0) = a;
            }
        }
    }
    else if (coord == 1)
    {
        const Real pi = 3.1415926535897932;
        if (dir == 0)
        {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real ri = offset[0] + dx[0]*(i);
                    Real a = std::abs((2.*pi)*ri*dx[1]);
                    area(i,j,0) = a;
                }
            }
        }
        else
        {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real rc = offset[0] + dx[0]*(i+0.5);
                    Real a = std::abs(dx[0]*(2.*pi)*rc);
                    area(i,j,0) = a;
                }
            }
        }
    }
    else
    {
        const Real pi = 3.1415926535897932;
        if (dir == 0)
        {
            for     (int j = lo.y; j <= hi.y; ++j) {
                Real ti = offset[1] + dx[1]*(j);
                Real to = ti + dx[1];
                Real tmp = (2.*pi)*(std::cos(ti)-std::cos(to));
                AMREX_PRAGMA_SIMD 
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real ri = offset[0] + dx[0]*(i);
                    Real a = tmp*ri*ri;
                    area(i,j,0) = a;
                }
            }
        }
        else
        {
            for     (int j = lo.y; j <= hi.y; ++j) {
                Real ti = offset[1] + dx[1]*(j);
                Real tmp = pi*std::sin(ti);
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real ri = offset[0] + dx[0]*(i);
                    Real ro = ri + dx[0];
                    Real a = tmp*(ro-ri)*(ro+ri);
                    area(i,j,0) = a;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
inline
void
amrex_setdloga (Box const& bx, Array4<Real> const& dloga,
                GpuArray<Real,2> const& offset,
                GpuArray<Real,2> const& dx, const int dir, const int coord) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    if (coord == 0)
    {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                dloga(i,j,0) = 0.0;
            }
        }
    }
    else if (coord == 1)
    {
        if (dir == 0)
        {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real rc = offset[0] + dx[0]*(i+0.5);
                    dloga(i,j,0) = 1.0/rc;
                }
            }
        }
        else
        {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    dloga(i,j,0) = 0.0;
                }
            }
        }
    }
    else
    {
        if (dir == 0)
        {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real rc = offset[0] + dx[0]*(i+0.5);
                    dloga(i,j,0) = 2.0/rc;
                }
            }
        }
        else
        {
            for     (int j = lo.y; j <= hi.y; ++j) {
                Real ti = offset[1] + dx[1]*(j);
                Real to = ti + dx[1];
                Real tmp = 1.0/std::tan(0.5*(ti+to));
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real rc = offset[0] + dx[0]*(i+0.5);
                    dloga(i,j,0) = tmp/rc;
                }
            }
        }
    }
}

}

#endif
