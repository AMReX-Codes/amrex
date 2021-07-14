#ifndef AMREX_MLNODELAP_K_H_
#define AMREX_MLNODELAP_K_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_LO_BCTYPES.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBCellFlag.H>
#endif

namespace amrex {
namespace {

    struct GetNode {
        AMREX_GPU_DEVICE Dim3 operator() (Dim3 const& lo, Dim3 const& len, int& offset)
        {
            Dim3 node;
            constexpr int nsten = AMREX_D_TERM(3,*3,*3);
            int icell = offset / nsten;
            node.z =  icell /        (len.x*len.y);
            node.y = (icell - node.z*(len.x*len.y)) /        len.x;
            node.x = (icell - node.z*(len.x*len.y)) - node.y*len.x;
            node.x += lo.x;
            node.y += lo.y;
            node.z += lo.z;
            offset -= icell*nsten;
            return node;
        }
    };

    struct GetNode2 {
        AMREX_GPU_DEVICE Dim3 operator() (int offset, Dim3 const& node)
        {
            // In 2D the offsets are
            //   6 7 8
            //   4 0 5
            //   1 2 3
            constexpr int nstenhalf = AMREX_SPACEDIM == 2 ? 4 : 13;
            if (offset == 0) {
                return node;
            } else {
                if (offset <= nstenhalf) { --offset; }
                Dim3 node2;
                node2.z = offset / 9;
                node2.y = (offset - node2.z*9) / 3;
                node2.x = (offset - node2.z*9) - node2.y*3;
                AMREX_D_TERM(node2.x += node.x-1;,
                             node2.y += node.y-1;,
                             node2.z += node.z-1);
                return node2;
            }
        }
    };

    constexpr int crse_cell = 0; // Do NOT change the values
    constexpr int fine_cell = 1;
    constexpr int crse_node = 0;
    constexpr int crse_fine_node = 1;
    constexpr int fine_node = 2;
#if (BL_USE_FLOAT)
    constexpr float eps = 1.e-30_rt;
#else
    constexpr double eps = 1.e-100_rt;
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
        if (lobc[idim] == LinOpBCType::Neumann || lobc[idim] == LinOpBCType::inflow) {
            Box const& blo = amrex::bdryLo(bx, idim);
            if (blo.smallEnd(idim) == nddom.smallEnd(idim)) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D (blo, i, j, k,
                {
                    rhs(i,j,k) *= s;
                });
            }
        }
        if (hibc[idim] == LinOpBCType::Neumann || hibc[idim] == LinOpBCType::inflow) {
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
    GpuArray<bool,AMREX_SPACEDIM> bflo{{AMREX_D_DECL(bclo[0] == LinOpBCType::Neumann ||
                                                     bclo[0] == LinOpBCType::inflow,
                                                     bclo[1] == LinOpBCType::Neumann ||
                                                     bclo[1] == LinOpBCType::inflow,
                                                     bclo[2] == LinOpBCType::Neumann ||
                                                     bclo[2] == LinOpBCType::inflow)}};
    GpuArray<bool,AMREX_SPACEDIM> bfhi{{AMREX_D_DECL(bchi[0] == LinOpBCType::Neumann ||
                                                     bchi[0] == LinOpBCType::inflow,
                                                     bchi[1] == LinOpBCType::Neumann ||
                                                     bchi[1] == LinOpBCType::inflow,
                                                     bchi[2] == LinOpBCType::Neumann ||
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
        if (!msk(i,j,k) && amrex::Math::abs(sten(i,j,k,0)) > s0_norm0) {
            x(i,j,k) /= sten(i,j,k,0);
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_sten (int i, int j, int k, Array4<Real> const& sol,
                          Real Ax, Array4<Real const> const& rhs,
                          Array4<Real const> const& sten,
                          Array4<int const> const& msk) noexcept
{
    if (msk(i,j,k)) {
        sol(i,j,k) = Real(0.0);
    } else if (sten(i,j,k,0) != Real(0.0)) {
        sol(i,j,k) += Real(2./3.) * (rhs(i,j,k) - Ax) / sten(i,j,k,0);
    }
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
            sol(i,j,k) = Real(0.0);
        } else if (sten(i,j,k,0) != Real(0.0)) {
            sol(i,j,k) += Real(2./3.) * (rhs(i,j,k) - Ax(i,j,k)) / sten(i,j,k,0);
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
