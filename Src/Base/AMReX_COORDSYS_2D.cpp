
#include <AMReX_COORDSYS_C.H>
#include <cmath>

namespace amrex {

AMREX_GPU_HOST_DEVICE
void amrex_setvol (Box const& bx, FArrayBox& vol,
                   GpuArray<Real,2> const& offset,
                   GpuArray<Real,2> const& dx, const int coord)
{
    const auto dp0 = vol.stridedPtr(bx);
    const IntVect len = bx.length();
    const IntVect lo  = bx.smallEnd();
    if (coord == 0) // Cartesian
    {
        Real dv = dx[0]*dx[1];
        for     (int j = 0; j < len[1]; ++j) {
            Real* AMREX_RESTRICT dp = dp0(0,j,0);
            for (int i = 0; i < len[0]; ++i) {
                dp[i] = dv;
            }
        }
    }
    else if (coord == 1) // r-z
    {
        const Real pi = 3.1415926535897932;
        for     (int j = 0; j < len[1]; ++j) {
            Real* AMREX_RESTRICT dp = dp0(0,j,0);
            for (int i = 0; i < len[0]; ++i) {
                Real ri = offset[0] + dx[0]*(i+lo[0]);
                Real ro = ri + dx[0];
                Real v = pi*dx[1]*dx[0]*(ro + ri);
                dp[i] = std::abs(v);
            }
        }
    }
    else // r-theta
    {
        const Real pi = 3.1415926535897932;
        for     (int j = 0; j < len[1]; ++j) {
            Real ti = offset[1] + dx[1]*(j+lo[1]);
            Real to = ti + dx[1];
            Real tmp = (2.*pi)*(std::cos(ti)-std::cos(to))/3.0;
            Real* AMREX_RESTRICT dp = dp0(0,j,0);
            for (int i = 0; i < len[0]; ++i) {
                Real ri = offset[0] + dx[0]*(i+lo[0]);
                Real ro = ri + dx[0];
                Real v = tmp*(ro-ri)*(ro*ro+ro*ri+ri*ri);
                dp[i] = std::abs(v);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_setarea (Box const& bx, FArrayBox& area,
                    GpuArray<Real,2> const& offset,
                    GpuArray<Real,2> const& dx, const int dir, const int coord)
{
    const auto dp0 = area.stridedPtr(bx);
    const IntVect len = bx.length();
    const IntVect lo  = bx.smallEnd();
    if (coord == 0)
    {
        Real a = (dir == 0) ? dx[1] : dx[0];
        for     (int j = 0; j < len[1]; ++j) {
            Real* AMREX_RESTRICT dp = dp0(0,j,0);
            for (int i = 0; i < len[0]; ++i) {
                dp[i] = a;
            }
        }
    }
    else if (coord == 1)
    {
        const Real pi = 3.1415926535897932;
        if (dir == 0)
        {
            for     (int j = 0; j < len[1]; ++j) {
                Real* AMREX_RESTRICT dp = dp0(0,j,0);
                for (int i = 0; i < len[0]; ++i) {
                    Real ri = offset[0] + dx[0]*(i+lo[0]);
                    Real a = std::abs((2.*pi)*ri*dx[1]);
                    dp[i] = a;
                }
            }
        }
        else
        {
            for     (int j = 0; j < len[1]; ++j) {
                Real* AMREX_RESTRICT dp = dp0(0,j,0);
                for (int i = 0; i < len[0]; ++i) {
                    Real rc = offset[0] + dx[0]*(i+lo[0]+0.5);
                    Real a = std::abs(dx[0]*(2.*pi)*rc);
                    dp[i] = a;
                }
            }
        }
    }
    else
    {
        const Real pi = 3.1415926535897932;
        if (dir == 0)
        {
            for     (int j = 0; j < len[1]; ++j) {
                Real ti = offset[1] + dx[1]*(j+lo[1]);
                Real to = ti + dx[1];
                Real tmp = (2.*pi)*(std::cos(ti)-std::cos(to));
                Real* AMREX_RESTRICT dp = dp0(0,j,0);
                for (int i = 0; i < len[0]; ++i) {
                    Real ri = offset[0] + dx[0]*(i+lo[0]);
                    Real a = tmp*ri*ri;
                    dp[i] = a;
                }
            }
        }
        else
        {
            for     (int j = 0; j < len[1]; ++j) {
                Real ti = offset[1] + dx[1]*(j+lo[1]);
                Real tmp = pi*std::sin(ti);
                Real* AMREX_RESTRICT dp = dp0(0,j,0);
                for (int i = 0; i < len[0]; ++i) {
                    Real ri = offset[0] + dx[0]*(i+lo[0]);
                    Real ro = ri + dx[0];
                    Real a = tmp*(ro-ri)*(ro+ri);
                    dp[i] = a;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_setdloga (Box const& bx, FArrayBox& dloga,
                     GpuArray<Real,2> const& offset,
                     GpuArray<Real,2> const& dx, const int dir, const int coord)
{
    const auto dp0 = dloga.stridedPtr(bx);
    const IntVect len = bx.length();
    const IntVect lo  = bx.smallEnd();
    if (coord == 0)
    {
        for     (int j = 0; j < len[1]; ++j) {
            Real* AMREX_RESTRICT dp = dp0(0,j,0);
            for (int i = 0; i < len[0]; ++i) {
                dp[i] = 0.0;
            }
        }
    }
    else if (coord == 1)
    {
        if (dir == 0)
        {
            for     (int j = 0; j < len[1]; ++j) {
                Real* AMREX_RESTRICT dp = dp0(0,j,0);
                for (int i = 0; i < len[0]; ++i) {
                    Real rc = offset[0] + dx[0]*(i+lo[0]+0.5);
                    dp[i] = 1.0/rc;
                }
            }
        }
        else
        {
            for     (int j = 0; j < len[1]; ++j) {
                Real* AMREX_RESTRICT dp = dp0(0,j,0);
                for (int i = 0; i < len[0]; ++i) {
                    dp[i] = 0.0;
                }
            }
        }
    }
    else
    {
        if (dir == 0)
        {
            for     (int j = 0; j < len[1]; ++j) {
                Real* AMREX_RESTRICT dp = dp0(0,j,0);
                for (int i = 0; i < len[0]; ++i) {
                    Real rc = offset[0] + dx[0]*(i+lo[0]+0.5);
                    dp[i] = 2.0/rc;
                }
            }
        }
        else
        {
            for     (int j = 0; j < len[1]; ++j) {
                Real ti = offset[1] + dx[1]*(j+lo[1]);
                Real to = ti + dx[1];
                Real tmp = 1.0/std::tan(0.5*(ti+to));
                Real* AMREX_RESTRICT dp = dp0(0,j,0);
                for (int i = 0; i < len[0]; ++i) {
                    Real rc = offset[0] + dx[0]*(i+lo[0]+0.5);
                    dp[i] = tmp/rc;
                }
            }
        }
    }
}

}
