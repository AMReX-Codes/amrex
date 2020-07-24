#ifndef AMREX_MLNODETENSORLAP_1D_K_H_
#define AMREX_MLNODETENSORLAP_1D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_interpadd (int i, int, int, Array4<Real> const& fine,
                          Array4<Real const> const& crse, Array4<int const> const& msk) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_adotx (Box const& b, Array4<Real> const& y, Array4<Real const> const& x,
                      GpuArray<Real,6> const& s,
                      GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_gauss_seidel (Box const& b, Array4<Real> const& sol,
                             Array4<Real const> const& rhs, Array4<int const> const& msk,
                             GpuArray<Real,6> const& s,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_normalize (Box const& b, Array4<Real> const& phi,
                          Array4<int const> const& msk, GpuArray<Real,6> const& s,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

}

#endif
