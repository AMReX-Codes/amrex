
#include <AMReX_COORDSYS_C.H>
#include <cmath>

namespace amrex {

AMREX_GPU_HOST_DEVICE
void amrex_setvol (Box const& bx, FArrayBox& vol,
                   GpuArray<Real,1> const& offset,
                   GpuArray<Real,1> const& dx, const int coord)
{
    const auto dp0 = vol.stridedPtr(bx);
    const IntVect len = bx.length();
    const int lo = bx.smallEnd(0);
    if (coord == 0) // Cartesian
    {
        Real dv = dx[0];
        Real* AMREX_RESTRICT dp = dp0(0,0,0);
        for (int i = 0; i < len[0]; ++i) {
            dp[i] = dv;
        }
    }
    else if (coord == 1)
    {
        const Real pi = 3.1415926535897932;
        Real* AMREX_RESTRICT dp = dp0(0,0,0);
        for (int i = 0; i < len[0]; ++i) {
            Real ri = offset[0] + dx[0]*(i+lo);
            Real ro = ri + dx[0];
            Real v = pi*(ro-ri)*(ro + ri);
            dp[i] = std::abs(v);
        }
    }
    else
    {
        const Real pi = 3.1415926535897932;
        Real* AMREX_RESTRICT dp = dp0(0,0,0);
        for (int i = 0; i < len[0]; ++i) {
            Real ri = offset[0] + dx[0]*(i+lo);
            Real ro = ri + dx[0];
            Real v = ((4./3.)*pi)*(ro-ri)*(ro*ro+ro*ri+ri*ri);
            dp[i] = std::abs(v);
        }
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_setarea (Box const& bx, FArrayBox& area,
                    GpuArray<Real,1> const& offset,
                    GpuArray<Real,1> const& dx, const int /*dir*/, const int coord)
{
    const auto dp0 = area.stridedPtr(bx);
    const IntVect len = bx.length();
    const int lo = bx.smallEnd(0);
    if (coord == 0)
    {
        Real* AMREX_RESTRICT dp = dp0(0,0,0);
        for (int i = 0; i < len[0]; ++i) {
            dp[i] = 1.0;
        }
    }
    else if (coord == 1)
    {
        const Real pi = 3.1415926535897932;
        Real* AMREX_RESTRICT dp = dp0(0,0,0);
        for (int i = 0; i < len[0]; ++i) {
            Real ri = offset[0] + dx[0]*(i+lo);
            Real a = (2.*pi)*ri;
            dp[i] = std::abs(a);
        }
    }
    else
    {
        const Real pi = 3.1415926535897932;
        Real* AMREX_RESTRICT dp = dp0(0,0,0);
        for (int i = 0; i < len[0]; ++i) {
            Real ri = offset[0] + dx[0]*(i+lo);
            Real a = (4.0*pi)*ri*ri;
            dp[i] = std::abs(a);
        }
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_setdloga (Box const& bx, FArrayBox& dloga,
                     GpuArray<Real,1> const& offset,
                     GpuArray<Real,1> const& dx, const int /*dir*/, const int coord)
{
    const auto dp0 = dloga.stridedPtr(bx);
    const IntVect len = bx.length();
    const int lo = bx.smallEnd(0);
    if (coord == 0)
    {
        Real* AMREX_RESTRICT dp = dp0(0,0,0);
        for (int i = 0; i < len[0]; ++i) {
            dp[i] = 0.0;
        }
    }
    else if (coord == 1)
    {
        Real* AMREX_RESTRICT dp = dp0(0,0,0);
        for (int i = 0; i < len[0]; ++i) {
            Real rc = offset[0] + dx[0]*(i+lo+0.5);
            dp[i] = 1.0/rc;
        }
    }
    else
    {
        Real* AMREX_RESTRICT dp = dp0(0,0,0);
        for (int i = 0; i < len[0]; ++i) {
            Real ri = offset[0] + dx[0]*(i+lo+0.5);
            dp[i] = 2.0/ri;
        }
    }
}

}
