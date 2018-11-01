
#include <AMReX_COORDSYS_C.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE
void amrex_setvol (Box const& bx, FArrayBox& vol,
                   GpuArray<Real,3> const& /*offset*/,
                   GpuArray<Real,3> const& dx, const int /*coord*/)
{
    Real dv = dx[0]*dx[1]*dx[2];
    const auto dp0 = vol.stridedPtr(bx);
    const IntVect len = bx.length();
    for         (int k = 0; k < len[2]; ++k) {
        for     (int j = 0; j < len[1]; ++j) {
            Real* AMREX_RESTRICT dp = dp0(0,j,k);
            for (int i = 0; i < len[0]; ++i) {
                dp[i] = dv;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_setarea (Box const& bx, FArrayBox& area,
                    GpuArray<Real,3> const& offset,
                    GpuArray<Real,3> const& dx, const int dir, const int /*coord*/)
{
    Real a = (dir == 0) ? dx[1]*dx[2] : ((dir == 1) ? dx[0]*dx[2] : dx[0]*dx[1]);
    const auto dp0 = area.stridedPtr(bx);
    const IntVect len = bx.length();
    for         (int k = 0; k < len[2]; ++k) {
        for     (int j = 0; j < len[1]; ++j) {
            Real* AMREX_RESTRICT dp = dp0(0,j,k);
            for (int i = 0; i < len[0]; ++i) {
                dp[i] = a;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_setdloga (Box const& bx, FArrayBox& dloga,
                     GpuArray<Real,3> const& /*offset*/,
                     GpuArray<Real,3> const& /*dx*/, const int /*dir*/, const int /*coord*/)
{
    const auto dp0 = dloga.stridedPtr(bx);
    const IntVect len = bx.length();
    for         (int k = 0; k < len[2]; ++k) {
        for     (int j = 0; j < len[1]; ++j) {
            Real* AMREX_RESTRICT dp = dp0(0,j,k);
            for (int i = 0; i < len[0]; ++i) {
                dp[i] = 0.0;
            }
        }
    }
}

}
