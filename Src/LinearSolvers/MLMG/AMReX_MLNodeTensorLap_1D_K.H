#ifndef AMREX_MLNODETENSORLAP_1D_K_H_
#define AMREX_MLNODETENSORLAP_1D_K_H_
#include <AMReX_Config.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_interpadd (int, int, int, Array4<Real> const&,
                          Array4<Real const> const&, Array4<int const> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_adotx (Box const&, Array4<Real> const&, Array4<Real const> const&,
                      GpuArray<Real,6> const&,
                      GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_gauss_seidel (Box const&, Array4<Real> const&,
                             Array4<Real const> const&, Array4<int const> const&,
                             GpuArray<Real,6> const&,
                             GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_normalize (Box const&, Array4<Real> const&,
                          Array4<int const> const&, GpuArray<Real,6> const&,
                          GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

}

#endif
