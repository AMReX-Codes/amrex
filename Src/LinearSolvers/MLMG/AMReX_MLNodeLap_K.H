#ifndef AMREX_MLNODELAP_K_H_
#define AMREX_MLNODELAP_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_LO_BCTYPES.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBCellFlag.H>
#endif

namespace amrex {
namespace {

    constexpr int crse_cell = 0;
    constexpr int fine_cell = 1;
    constexpr int crse_node = 0;
    constexpr int crse_fine_node = 1;
    constexpr int fine_node = 2;
#if (BL_USE_FLOAT)
    constexpr float eps = 1.e-30;
#else
    constexpr double eps = 1.e-100;
#endif
    constexpr Real almostone = Real(1.) - Real(100.)*std::numeric_limits<Real>::epsilon();

}

namespace {
inline void
mlndlap_scale_neumann_bc (Real s, Box const& bx, Array4<Real> const& rhs, Box const& nddom,
                          GpuArray<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                          GpuArray<LinOpBCType,AMREX_SPACEDIM> const& hibc) noexcept
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (lobc[idim] == LinOpBCType::Neumann or lobc[idim] == LinOpBCType::inflow) {
            Box const& blo = amrex::bdryLo(bx, idim);
            if (blo.smallEnd(idim) == nddom.smallEnd(idim)) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D (blo, i, j, k,
                {
                    rhs(i,j,k) *= s;
                });
            }
        }
        if (hibc[idim] == LinOpBCType::Neumann or hibc[idim] == LinOpBCType::inflow) {
            Box const& bhi = amrex::bdryHi(bx, idim);
            if (bhi.bigEnd(idim) == nddom.bigEnd(idim)) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bhi, i, j, k,
                {
                    rhs(i,j,k) *= s;
                });
            }
        }
    }
}
}

inline void
mlndlap_impose_neumann_bc (Box const& bx, Array4<Real> const& rhs, Box const& nddom,
                           GpuArray<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                           GpuArray<LinOpBCType,AMREX_SPACEDIM> const& hibc) noexcept
{
    mlndlap_scale_neumann_bc(2.0, bx, rhs, nddom, lobc, hibc);
}

inline void
mlndlap_unimpose_neumann_bc (Box const& bx, Array4<Real> const& rhs, Box const& nddom,
                           GpuArray<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                           GpuArray<LinOpBCType,AMREX_SPACEDIM> const& hibc) noexcept
{
    mlndlap_scale_neumann_bc(0.5, bx, rhs, nddom, lobc, hibc);
}

}

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLNodeLap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLNodeLap_2D_K.H>
#else
#include <AMReX_MLNodeLap_3D_K.H>
#endif

namespace amrex {

template <typename T>
void mlndlap_fillbc_cc (Box const& vbx, Array4<T> const& sigma, Box const& domain,
                        GpuArray<LinOpBCType, AMREX_SPACEDIM> bclo,
                        GpuArray<LinOpBCType, AMREX_SPACEDIM> bchi) noexcept
{
    GpuArray<bool,AMREX_SPACEDIM> bflo{{AMREX_D_DECL(bclo[0] != LinOpBCType::Periodic,
                                                     bclo[1] != LinOpBCType::Periodic,
                                                     bclo[2] != LinOpBCType::Periodic)}};
    GpuArray<bool,AMREX_SPACEDIM> bfhi{{AMREX_D_DECL(bchi[0] != LinOpBCType::Periodic,
                                                     bchi[1] != LinOpBCType::Periodic,
                                                     bchi[2] != LinOpBCType::Periodic)}};
    mlndlap_bc_doit(vbx, sigma, domain, bflo, bfhi);
}

template <typename T>
void mlndlap_applybc (Box const& vbx, Array4<T> const& phi, Box const& domain,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> bclo,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> bchi) noexcept
{
    GpuArray<bool,AMREX_SPACEDIM> bflo{{AMREX_D_DECL(bclo[0] == LinOpBCType::Neumann or
                                                     bclo[0] == LinOpBCType::inflow,
                                                     bclo[1] == LinOpBCType::Neumann or
                                                     bclo[1] == LinOpBCType::inflow,
                                                     bclo[2] == LinOpBCType::Neumann or
                                                     bclo[2] == LinOpBCType::inflow)}};
    GpuArray<bool,AMREX_SPACEDIM> bfhi{{AMREX_D_DECL(bchi[0] == LinOpBCType::Neumann or
                                                     bchi[0] == LinOpBCType::inflow,
                                                     bchi[1] == LinOpBCType::Neumann or
                                                     bchi[1] == LinOpBCType::inflow,
                                                     bchi[2] == LinOpBCType::Neumann or
                                                     bchi[2] == LinOpBCType::inflow)}};
    mlndlap_bc_doit(vbx, phi, domain, bflo, bfhi);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_sten (Box const& bx, Array4<Real> const& x,
                             Array4<Real const> const& sten,
                             Array4<int const> const& msk, Real s0_norm0) noexcept
{
    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (!msk(i,j,k) and amrex::Math::abs(sten(i,j,k,0)) > s0_norm0) {
            x(i,j,k) /= sten(i,j,k,0);
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_sten (Box const& bx, Array4<Real> const& sol,
                          Array4<Real const> const& Ax,
                          Array4<Real const> const& rhs,
                          Array4<Real const> const& sten,
                          Array4<int const> const& msk) noexcept
{
    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else if (sten(i,j,k,0) != 0.0) {
            sol(i,j,k) += (2./3.) * (rhs(i,j,k) - Ax(i,j,k)) / sten(i,j,k,0);
        }
    });
}

AMREX_FORCE_INLINE
bool mlndlap_any_fine_sync_cells (Box const& bx, Array4<int const> const& msk, int fine_flag) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        if (msk(i,j,k) == fine_flag) return true;
    }}}
    return false;
}

}

#endif
