#ifndef AMREX_MLNODELAP_1D_K_H_
#define AMREX_MLNODELAP_1D_K_H_
#include <AMReX_Config.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_nodal_mask (int /*i*/, int /*j*/, int /*k*/, Array4<int> const& /*nmsk*/,
                             Array4<int const> const& /*cmsk*/) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dirichlet_mask (Box const&, Array4<int> const&,
                                 Array4<int const> const&, Box const&,
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const&,
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dot_mask (Box const&, Array4<Real> const&,
                           Array4<int const> const&, Box const&,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const&,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_zero_fine (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                        Array4<int const> const&, int) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_avgdown_coeff_x (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                              Array4<Real const> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_avgdown_coeff (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                                 Array4<Real const> const&, int) noexcept
{}

template <typename T>
inline void mlndlap_bc_doit (Box const&, Array4<T> const&, Box const&,
                             GpuArray<bool,AMREX_SPACEDIM> const&,
                             GpuArray<bool,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndlap_adotx_c (int /*i*/, int /*j*/, int /*k*/, Array4<Real const> const&,
                      Real, Array4<int const> const&,
                      GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{ return Real(0.0); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndlap_adotx_ha (int /*i*/, int /*j*/, int /*k*/, Array4<Real const> const&,
                       Array4<Real const> const&, Array4<int const> const&,
                       GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{ return Real(0.0); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndlap_adotx_aa (int /*i*/, int /*j*/, int /*k*/, Array4<Real const> const&,
                       Array4<Real const> const&, Array4<int const> const&,
                       GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{ return Real(0.0); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_ha (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&, Array4<Real const> const&,
                           Array4<int const> const&, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_aa (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&, Array4<Real const> const&,
                           Array4<int const> const&, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_ha (int, int, int, Array4<Real> const&, Real,
                        Array4<Real const> const&, Array4<Real const> const&,
                        Array4<int const> const&, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_ha (Box const&, Array4<Real> const&, Array4<Real const> const&,
                        Array4<Real const> const&, Array4<Real const> const&,
                        Array4<int const> const&, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_aa (int, int, int, Array4<Real> const&, Real,
                        Array4<Real const> const&, Array4<Real const> const&,
                        Array4<int const> const&, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_aa (Box const&, Array4<Real> const&, Array4<Real const> const&,
                        Array4<Real const> const&, Array4<Real const> const&,
                        Array4<int const> const&, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_c (int, int, int, Array4<Real> const&, Real,
                       Array4<Real const> const&, Real,
                       Array4<int const> const&, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_c (Box const&, Array4<Real> const&, Array4<Real const> const&,
                       Array4<Real const> const&, Real,
                       Array4<int const> const&, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_ha (Box const&, Array4<Real> const&,
                              Array4<Real const> const&, Array4<Real const> const&,
                              Array4<int const> const&, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_aa (Box const&, Array4<Real> const&,
                              Array4<Real const> const&, Array4<Real const> const&,
                              Array4<int const> const&, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_c (Box const&, Array4<Real> const&,
                             Array4<Real const> const&, Real,
                             Array4<int const> const&, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_with_line_solve_aa(Box const&, Array4<Real> const&,
                                             Array4<Real const> const&, Array4<Real const> const&,
                                             Array4<int const> const&, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                          Array4<Real const> const&, Array4<int const> const&) noexcept
{}

template <int rr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                          Array4<Real const> const&, Array4<int const> const&,
                          Box const&,
                          GpuArray<LinOpBCType, AMREX_SPACEDIM> const&,
                          GpuArray<LinOpBCType, AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_restriction (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                          Array4<Real const> const&, Array4<int const> const&, int) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_c (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                          Array4<Real const> const&,
                          Array4<int const> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_aa (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                           Array4<Real const> const&, Array4<Real const> const&,
                           Array4<int const> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_interpadd_aa (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                           Array4<Real const> const&, Array4<Real const> const&,
                           Array4<int const> const&, int) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_ha (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                           Array4<Real const> const&, Array4<Real const> const&,
                           Array4<int const> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&, Array4<Real const> const&,
                   Array4<int const> const&,
                   GpuArray<Real,AMREX_SPACEDIM> const&,
                   Box const&,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const&,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndlap_rhcc (int /*i*/, int /*j*/, int /*k*/, Array4<Real const> const&,
                   Array4<int const> const&) noexcept
{ return Real(0.); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&, Array4<Real const> const&,
                     Array4<Real const> const&, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu_c (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&, Array4<Real const> const&,
                     Real, GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

template <int rr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_rhcc_fine_contrib (int /*i*/, int /*j*/, int /*k*/, Box const&,
                                Array4<Real> const&, Array4<Real const> const&,
                                Array4<int const> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_cf_contrib (int, int, int, Array4<Real> const&,
                              Array4<Real const> const&, Array4<Real const> const&,
                              Array4<Real const> const&, Array4<int const> const&,
                              Array4<int const> const&, Array4<int const> const&,
                              GpuArray<Real,AMREX_SPACEDIM> const&,
                              Box const&, GpuArray<LinOpBCType,AMREX_SPACEDIM> const&,
                              GpuArray<LinOpBCType,AMREX_SPACEDIM> const&,
                              bool) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_crse_resid (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                         Array4<Real const> const&, Array4<int const> const&,
                         Box const&, GpuArray<LinOpBCType,AMREX_SPACEDIM> const&,
                         GpuArray<LinOpBCType,AMREX_SPACEDIM> const&,
                         bool) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_cf_contrib (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                             Array4<Real const> const&, Array4<Real const> const&,
                             Array4<Real const> const&, Array4<int const> const&,
                             Array4<int const> const&, Array4<int const> const&,
                             Array4<Real const> const&,
                             GpuArray<Real,AMREX_SPACEDIM> const&, Box const&,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const&,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const&,
                             bool) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_cf_contrib_cs (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                                Array4<Real const> const&, Array4<Real const> const&,
                                Real, Array4<int const> const&,
                                Array4<int const> const&, Array4<int const> const&,
                                Array4<Real const> const&,
                                GpuArray<Real,AMREX_SPACEDIM> const&, Box const&,
                                GpuArray<LinOpBCType, AMREX_SPACEDIM> const&,
                                GpuArray<LinOpBCType, AMREX_SPACEDIM> const&,
                                bool) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil (Box const&, Array4<Real> const&,
                          Array4<Real const> const&,
                          GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil_s0 (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_stencil_rap (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                          Array4<Real const> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndlap_adotx_sten (int /*i*/, int /*j*/, int /*k*/, Array4<Real const> const&,
                         Array4<Real const> const&, Array4<int const> const&) noexcept
{ return Real(0.0); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_sten (Box const&, Array4<Real> const&,
                                Array4<Real const> const&,
                                Array4<Real const> const&,
                                Array4<int const> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_rap (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                            Array4<Real const> const&, Array4<Real const> const&,
                            Array4<int const> const&) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction_rap (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                              Array4<Real const> const&, Array4<Real const> const&,
                              Array4<int const> const&) noexcept
{}

}

#endif
