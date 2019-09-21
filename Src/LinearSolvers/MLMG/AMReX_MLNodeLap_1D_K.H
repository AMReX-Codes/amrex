#ifndef AMREX_MLNODELAP_1D_K_H_
#define AMREX_MLNODELAP_1D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu (int i, int j, int k, Array4<Real> const& rhs, Array4<Real const> const& vel,
                   Array4<int const> const& msk,
                   GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

}
#endif
