#ifndef _Adv_3d_cpp_
#define _Adv_3d_cpp_

#include <Adv_3d.H>
#include <AmrCoreAdv_F.H>
#include <AMReX_FArrayBox.H>

using namespace amrex;

AMREX_GPU_DEVICE
void conservative(Box const& bx,
                  const FArrayBox& statein,
                  FArrayBox& stateout,
                  AMREX_D_DECL(FArrayBox& fx,
                               FArrayBox& fy,
                               FArrayBox& fz),
                  const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);

    const auto uin  = statein.view(lo);
    const auto uout = stateout.view(lo);
    const auto flxx = fx.view(lo);
    const auto flxy = fy.view(lo);
    const auto flxz = fz.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                uout(i,j,k) = uin(i,j,k) + 
                    ( (flxx(i,j,k) - flxx(i+1,j,k)) * dtdx[0] 
                    + (flxy(i,j,k) - flxy(i,j+1,k)) * dtdx[1] 
                    + (flxz(i,j,k) - flxz(i,j,k+1)) * dtdx[2] );
            }
        }
    }
}

AMREX_GPU_DEVICE
void flux_scale_x(Box const& bx,
                    FArrayBox& fx,
                    const Real& dt,
                    const GpuArray<Real, AMREX_SPACEDIM>& dx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto flxx = fx.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                flxx(i,j,k) = flxx(i,j,k) * (dt * dx[1]*dx[2]);
            }
        }
    }
}

AMREX_GPU_DEVICE
void flux_scale_y(Box const& bx,
                    FArrayBox& fy,
                    const Real& dt,
                    const GpuArray<Real, AMREX_SPACEDIM>& dx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto flxy = fy.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                flxy(i,j,k) = flxy(i,j,k) * (dt * dx[0]*dx[2]);
            }
        }
    }
}

AMREX_GPU_DEVICE
void flux_scale_z(Box const& bx,
                    FArrayBox& fz,
                    const Real& dt,
                    const GpuArray<Real, AMREX_SPACEDIM>& dx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto flxz = fz.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                flxz(i,j,k) = flxz(i,j,k) * (dt * dx[0]*dx[1]);
            }
        }
    }
}

#endif
