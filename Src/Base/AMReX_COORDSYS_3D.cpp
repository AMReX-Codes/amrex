
#include <AMReX_COORDSYS_C.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE
void amrex_setvol (Box const& bx, FArrayBox& vol,
                   GpuArray<Real,3> const& /*offset*/,
                   GpuArray<Real,3> const& dx, const int /*coord*/)
{
    Real dv = dx[0]*dx[1]*dx[2];
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = vol.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                dp(i,j,k) = dv;
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

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = area.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                dp(i,j,k) = a;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE
void amrex_setdloga (Box const& bx, FArrayBox& dloga,
                     GpuArray<Real,3> const& /*offset*/,
                     GpuArray<Real,3> const& /*dx*/, const int /*dir*/, const int /*coord*/)
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = dloga.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                dp(i,j,k) = 0.0;
            }
        }
    }
}

}
