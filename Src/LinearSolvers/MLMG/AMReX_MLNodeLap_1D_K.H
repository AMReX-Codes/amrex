#ifndef AMREX_MLNODELAP_1D_K_H_
#define AMREX_MLNODELAP_1D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_nodal_mask (int /*i*/, int /*j*/, int /*k*/, Array4<int> const& /*nmsk*/,
                             Array4<int const> const& /*cmsk*/) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dirichlet_mask (Box const& bx, Array4<int> const& dmsk,
                                 Array4<int const> const& omsk, Box const& dom,
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dot_mask (Box const& bx, Array4<Real> const& dmsk,
                           Array4<int const> const& omsk, Box const& dom,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_zero_fine (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& phi,
                        Array4<int const> const& msk, int fine_flag) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_avgdown_coeff_x (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& crse,
                              Array4<Real const> const& fine) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_avgdown_coeff (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& crse,
                              Array4<Real const> const& fine, int idir) noexcept
{}

template <typename T>
inline void mlndlap_bc_doit (Box const& vbx, Array4<T> const& a, Box const& domain,
                             GpuArray<bool,AMREX_SPACEDIM> const& bflo,
                             GpuArray<bool,AMREX_SPACEDIM> const& bfhi) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_ha (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& y, Array4<Real const> const& x,
                       Array4<Real const> const& sx, Array4<int const> const& msk,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_aa (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& y, Array4<Real const> const& x,
                       Array4<Real const> const& sig, Array4<int const> const& msk,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_ha (Box const& bx, Array4<Real> const& x, Array4<Real const> const& sx,
                           Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_aa (Box const& bx, Array4<Real> const& x, Array4<Real const> const& sig,
                           Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_ha (Box const& bx, Array4<Real> const& sol, Array4<Real const> const& Ax,
                        Array4<Real const> const& rhs, Array4<Real const> const& sx,
                        Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_aa (Box const& bx, Array4<Real> const& sol, Array4<Real const> const& Ax,
                        Array4<Real const> const& rhs, Array4<Real const> const& sig,
                        Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_ha (Box const& bx, Array4<Real> const& sol,
                              Array4<Real const> const& rhs, Array4<Real const> const& sx,
                              Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_aa (Box const& bx, Array4<Real> const& sol,
                              Array4<Real const> const& rhs, Array4<Real const> const& sig,
                              Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_with_line_solve_aa(Box const& bx, Array4<Real> const& sol,
                                             Array4<Real const> const& rhs, Array4<Real const> const& sig,
                                             Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& crse,
                          Array4<Real const> const& fine, Array4<int const> const& msk) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_restriction (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& crse,
                          Array4<Real const> const& fine, Array4<int const> const& msk, int idir) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_aa (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& fine,
                           Array4<Real const> const& crse, Array4<Real const> const& sig,
                           Array4<int const> const& msk) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_interpadd_aa (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& fine,
                           Array4<Real const> const& crse, Array4<Real const> const& sig,
                           Array4<int const> const& msk, int idir) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_ha (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& fine,
                           Array4<Real const> const& crse, Array4<Real const> const& sigx,
                           Array4<int const> const& msk) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& rhs, Array4<Real const> const& vel,
                   Array4<int const> const& msk,
                   GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                   Box const& nddom,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndlap_rhcc (int /*i*/, int /*j*/, int /*k*/, Array4<Real const> const& rhcc,
                   Array4<int const> const& msk) noexcept
{ return Real(0.); }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& u, Array4<Real const> const& p,
                     Array4<Real const> const& sig, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_compute_fine_contrib (int /*i*/, int /*j*/, int /*k*/, Box const& fvbx,
                                        Array4<Real> const& frh, Array4<Real const> const& vel,
                                        GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                        Box const& nddom,
                                        GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                        GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_add_fine_contrib (int /*i*/, int /*j*/, int /*k*/, Box const& fvbx,
                                    Array4<Real> const& rhs, Array4<Real const> const& frh,
                                    Array4<int const> const& msk) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_rhcc_fine_contrib (int /*i*/, int /*j*/, int /*k*/, Box const& fvbx,
                                Array4<Real> const& rhs, Array4<Real const> const& cc,
                                Array4<int const> const& msk) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_cf_contrib (int i, int j, int, Array4<Real> const& rhs,
                              Array4<Real const> const& vel, Array4<Real const> const& fc,
                              Array4<Real const> const& rhcc, Array4<int const> const& dmsk,
                              Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Box const& nddom, GpuArray<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                              GpuArray<LinOpBCType,AMREX_SPACEDIM> const& hibc,
                              bool neumann_doubling) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_crse_resid (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& resid,
                         Array4<Real const> const& rhs, Array4<int const> const& msk,
                         Box const& nddom, GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                         GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi,
                         bool neumann_doubling) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_fine_Ax (int /*i*/, int /*j*/, int /*k*/, Box const& fvbx, Array4<Real> const& Ax,
                          Array4<Real const> const& x, Array4<Real const> const& sig,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_fine_contrib (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& f,
                               Array4<Real const> const& Ax, Array4<int const> const& msk) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_cf_contrib (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& res,
                             Array4<Real const> const& phi, Array4<Real const> const& rhs,
                             Array4<Real const> const& sig, Array4<int const> const& dmsk,
                             Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                             Array4<Real const> const& fc,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv, Box const& nddom,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& lobc,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& hibc,
                             bool neumann_doubling) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil (Box const& bx, Array4<Real> const& sten,
                          Array4<Real const> const& sig,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil_s0 (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& sten) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_stencil_rap (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& csten,
                          Array4<Real const> const& fsten) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_sten (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& y, Array4<Real const> const& x,
                         Array4<Real const> const& sten, Array4<int const> const& msk) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_sten (Box const& bx, Array4<Real> const& sol,
                                Array4<Real const> const& rhs,
                                Array4<Real const> const& sten,
                                Array4<int const> const& msk) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_rap (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& fine,
                            Array4<Real const> const& crse, Array4<Real const> const& sten,
                            Array4<int const> const& msk) noexcept
{}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction_rap (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& crse,
                              Array4<Real const> const& fine, Array4<Real const> const& sten,
                              Array4<int const> const& msk) noexcept
{}

}

#endif
