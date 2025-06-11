#ifndef AMREX_MLNODEABECLAP_1D_K_H_
#define AMREX_MLNODEABECLAP_1D_K_H_

namespace amrex {

inline void
mlndabeclap_gauss_seidel_aa (Box const& /*bx*/, Array4<Real> const& /*sol*/,
                             Array4<Real const> const& /*rhs*/,
                             Real /*alpha*/, Real /*beta*/,
                             Array4<Real const> const& /*acf*/,
                             Array4<Real const> const& /*bcf*/,
                             Array4<int const> const& /*msk*/,
                             GpuArray<Real,AMREX_SPACEDIM> const& /*dxinv*/) noexcept
{}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
mlndabeclap_jacobi_aa (int /*i*/, int /*j*/, int /*k*/,
                       Array4<Real> const& /*sol*/,
                       Real /*lap*/,
                       Array4<Real const> const& /*rhs*/,
                       Real /*alpha*/, Real /*beta*/,
                       Array4<Real const> const& /*acf*/,
                       Array4<Real const> const& /*bcf*/,
                       Array4<int const> const& /*msk*/,
                       GpuArray<Real,AMREX_SPACEDIM> const& /*dxinv*/) noexcept
{}

}

#endif
