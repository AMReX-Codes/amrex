#ifndef Adv_3d_H_
#define Adv_3d_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

using namespace amrex;

AMREX_GPU_DEVICE
void conservative(Box const& bx,
                  const FArrayBox& statein,
                  FArrayBox& stateout,
                  AMREX_D_DECL(FArrayBox& fx,
                               FArrayBox& fy,
                               FArrayBox& fz),
                  const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

AMREX_GPU_DEVICE
void flux_scale_x(Box const& bx,
                  FArrayBox& fx,
                  const Real& dt,
                  const GpuArray<Real, AMREX_SPACEDIM>& dx);

AMREX_GPU_DEVICE
void flux_scale_y(Box const& bx,
                  FArrayBox& fy,
                  const Real& dt,
                  const GpuArray<Real, AMREX_SPACEDIM>& dx);

AMREX_GPU_DEVICE
void flux_scale_z(Box const& bx,
                  FArrayBox& fz,
                  const Real& dt,
                  const GpuArray<Real, AMREX_SPACEDIM>& dx);

#endif
