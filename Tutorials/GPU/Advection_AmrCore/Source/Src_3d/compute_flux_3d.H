#ifndef compute_flux_3d_H
#define compute_flux_3d_H

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

using namespace amrex;

AMREX_GPU_DEVICE
void flux_x(Box const& bx,
            const FArrayBox& state,
            const FArrayBox& velx,
            FArrayBox& phix,
            const FArrayBox& slope,
            const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

AMREX_GPU_DEVICE
void flux_y(Box const& bx,
            const FArrayBox& state,
            const FArrayBox& vely,
            FArrayBox& phiy,
            const FArrayBox& slope,
            const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

AMREX_GPU_DEVICE
void flux_z(Box const& bx,
            const FArrayBox& state,
            const FArrayBox& velz,
            FArrayBox& phiz,
            const FArrayBox& slope,
            const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

AMREX_GPU_DEVICE
void flux_xy(Box const& bx,
             AMREX_D_DECL(const FArrayBox& velx,
                          const FArrayBox& vely,
                          const FArrayBox& velz),
             AMREX_D_DECL(const FArrayBox& phix,
                          const FArrayBox& phiy,
                          const FArrayBox& phiz),
             FArrayBox& phix_y,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

AMREX_GPU_DEVICE
void flux_xz(Box const& bx,
             AMREX_D_DECL(const FArrayBox& velx,
                          const FArrayBox& vely,
                          const FArrayBox& velz),
             AMREX_D_DECL(const FArrayBox& phix,
                          const FArrayBox& phiy,
                          const FArrayBox& phiz),
             FArrayBox& phix_z,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

AMREX_GPU_DEVICE
void flux_yx(Box const& bx,
             AMREX_D_DECL(const FArrayBox& velx,
                          const FArrayBox& vely,
                          const FArrayBox& velz),
             AMREX_D_DECL(const FArrayBox& phix,
                          const FArrayBox& phiy,
                          const FArrayBox& phiz),
             FArrayBox& phiy_x,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

AMREX_GPU_DEVICE
void flux_yz(Box const& bx,
             AMREX_D_DECL(const FArrayBox& velx,
                          const FArrayBox& vely,
                          const FArrayBox& velz),
             AMREX_D_DECL(const FArrayBox& phix,
                          const FArrayBox& phiy,
                          const FArrayBox& phiz),
             FArrayBox& phi_yz,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

AMREX_GPU_DEVICE
void flux_zx(Box const& bx,
             AMREX_D_DECL(const FArrayBox& velx,
                          const FArrayBox& vely,
                          const FArrayBox& velz),
             AMREX_D_DECL(const FArrayBox& phix,
                          const FArrayBox& phiy,
                          const FArrayBox& phiz),
             FArrayBox& phi_zx,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

AMREX_GPU_DEVICE
void flux_zy(Box const& bx,
             AMREX_D_DECL(const FArrayBox& velx,
                          const FArrayBox& vely,
                          const FArrayBox& velz),
             AMREX_D_DECL(const FArrayBox& phix,
                          const FArrayBox& phiy,
                          const FArrayBox& phiz),
             FArrayBox& phi_zy,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

AMREX_GPU_DEVICE
void combine_flux_x(Box const& bx,
                    const FArrayBox& velx,
                    const FArrayBox& vely,
                    const FArrayBox& velz,
                    FArrayBox& phix,
                    const FArrayBox& phiy_z,
                    const FArrayBox& phiz_y,
                    FArrayBox& flxx,
                    const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

AMREX_GPU_DEVICE
void combine_flux_y(Box const& bx,
                    const FArrayBox& velx,
                    const FArrayBox& vely,
                    const FArrayBox& velz,
                    FArrayBox& phiy,
                    const FArrayBox& phix_z,
                    const FArrayBox& phiz_x,
                    FArrayBox& flxy,
                    const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

AMREX_GPU_DEVICE
void combine_flux_z(Box const& bx,
                    const FArrayBox& velx,
                    const FArrayBox& vely,
                    const FArrayBox& velz,
                    FArrayBox& phiz,
                    const FArrayBox& phix_y,
                    const FArrayBox& phiy_x,
                    FArrayBox& flxz,
                    const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

#endif
