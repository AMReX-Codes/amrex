#ifndef AMREX_MLNODELAP_K_H_
#define AMREX_MLNODELAP_K_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_LO_BCTYPES.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBCellFlag.H>
#endif
#include <AMReX_MLNodeLinOp_K.H>

#include <limits>
#include <type_traits>

namespace amrex {

/// \cond DOXYGEN_IGNORE
namespace nodelap_detail {

    template <typename T>
    inline constexpr T eps_v = T(1.e-100L);
    template <>
    inline constexpr float eps_v<float> = float(1.e-30L);

    template <typename T>
    inline constexpr T almostone_v = T(1.) - T(100.)*std::numeric_limits<T>::epsilon();
    template <typename T>
    inline constexpr T almostzero_v = T(1.) - almostone_v<T>;

    constexpr Real eps = eps_v<Real>;
    constexpr Real almostone  = almostone_v<Real>;
    constexpr Real almostzero = almostzero_v<Real>;

}
/// \endcond

template <typename T>
inline void
mlndlap_scale_neumann_bc (std::type_identity_t<T> s, Box const& bx, Array4<T> const& rhs, Box const& nddom,
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

template <typename T>
inline void
mlndlap_impose_neumann_bc (Box const& bx, Array4<T> const& rhs, Box const& nddom,
                           GpuArray<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                           GpuArray<LinOpBCType,AMREX_SPACEDIM> const& hibc) noexcept
{
    mlndlap_scale_neumann_bc(2.0, bx, rhs, nddom, lobc, hibc);
}

template <typename T>
inline void
mlndlap_unimpose_neumann_bc (Box const& bx, Array4<T> const& rhs, Box const& nddom,
                           GpuArray<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                           GpuArray<LinOpBCType,AMREX_SPACEDIM> const& hibc) noexcept
{
    mlndlap_scale_neumann_bc(0.5, bx, rhs, nddom, lobc, hibc);
}

namespace nodelap_detail {

// The Neumann (wall) domain faces a node lies on: wall[d] = 0 (lo), 1 (hi) or -1 (none).
// Returns the number of walls.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int mlndlap_walls (IntVect const& iv, Box const& nddom,
                   GpuArray<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                   GpuArray<LinOpBCType,AMREX_SPACEDIM> const& hibc,
                   GpuArray<int,AMREX_SPACEDIM>& wall) noexcept
{
    int nwall = 0;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        wall[d] = -1;
        if (lobc[d] == LinOpBCType::Neumann && iv[d] == nddom.smallEnd(d)) { wall[d] = 0; ++nwall; }
        if (hibc[d] == LinOpBCType::Neumann && iv[d] == nddom.bigEnd(d))   { wall[d] = 1; ++nwall; }
    }
    return nwall;
}

// Cell index along direction d for the node iv: row m (0, 1) from the wall if d is a wall
// direction, else the cell on side s (0: low, 1: high) of the node.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int mlndlap_cell_index (int ivd, int walld, int m, int s) noexcept
{
    if (walld == 0) { return ivd + m; }
    if (walld == 1) { return ivd - 1 - m; }
    return ivd - 1 + s;
}

// The stencil's contribution of velocity component d at a wall node (after the Neumann
// scaling), with wall direction t taken at row mt and every other wall direction at row 0:
// for d normal to a wall, 2/h times u_d in the first cell; otherwise the derivative
// (u_d(hi) - u_d(lo))/h across the node; both averaged over the cells in the non-wall
// directions.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndlap_wall_term (IntVect const& iv, int d, int t, int mt,
                        GpuArray<int,AMREX_SPACEDIM> const& wall,
                        Array4<Real const> const& vel,
                        GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real sum = Real(0.0);
    int cnt = 0;
    // Enumerate the cells: for non-wall directions other than d, both sides.
    int ns[AMREX_SPACEDIM];
    for (int o = 0; o < AMREX_SPACEDIM; ++o) {
        ns[o] = (wall[o] < 0 && o != d) ? 2 : 1;
    }
    int const ncomb = AMREX_D_TERM(ns[0], *ns[1], *ns[2]);
    for (int c = 0; c < ncomb; ++c) {
        int rem = c;
        IntVect cell;
        for (int o = 0; o < AMREX_SPACEDIM; ++o) {
            int const sidx = rem % ns[o];
            rem /= ns[o];
            int const m = (o == t) ? mt : 0;
            cell[o] = mlndlap_cell_index(iv[o], wall[o], m, sidx);
        }
        if (wall[d] >= 0) {
            // Flux out of the control volume through its face h/2 from the wall: + at a
            // lo wall, - at a hi wall, as in mlndlap_divu.
            Real const sgn = (wall[d] == 0) ? Real(1.0) : Real(-1.0);
            sum += sgn * Real(2.0) * dxinv[d] * vel(cell, d);
        } else {
            IntVect clo(cell), chi(cell);
            clo[d] = iv[d] - 1;
            chi[d] = iv[d];
            sum += (vel(chi, d) - vel(clo, d)) * dxinv[d];
        }
        ++cnt;
    }
    return sum / Real(cnt);
}

// The stencil's source term at a wall node (after the Neumann scaling): the mean of rhcc
// over the cells touching the node, with wall direction t at row mt and the other wall
// directions at row 0.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndlap_wall_src (IntVect const& iv, int t, int mt,
                       GpuArray<int,AMREX_SPACEDIM> const& wall,
                       Array4<Real const> const& rhcc) noexcept
{
    Real sum = Real(0.0);
    int cnt = 0;
    int ns[AMREX_SPACEDIM];
    for (int o = 0; o < AMREX_SPACEDIM; ++o) {
        ns[o] = (wall[o] < 0) ? 2 : 1;
    }
    int const ncomb = AMREX_D_TERM(ns[0], *ns[1], *ns[2]);
    for (int c = 0; c < ncomb; ++c) {
        int rem = c;
        IntVect cell;
        for (int o = 0; o < AMREX_SPACEDIM; ++o) {
            int const sidx = rem % ns[o];
            rem /= ns[o];
            int const m = (o == t) ? mt : 0;
            cell[o] = mlndlap_cell_index(iv[o], wall[o], m, sidx);
        }
        sum += rhcc(cell);
        ++cnt;
    }
    return sum / Real(cnt);
}

}

// Second-order nodal divergence at walls.  The control volume of a wall node is half a
// cell in each wall direction (a quarter cell at a 2D corner).  mlndlap_divu evaluates
// the fluxes through its faces, and the tangential derivatives, at the first cell centres,
// h/2 from each wall; their mean over the control volume is their value at h/4 to second
// order, which differs by (h/4) d/dn.  This adds, after the Neumann scaling, the linear
// extrapolation from the first two rows of cells in every wall direction t that is not the
// component's own direction: 1/4 [T(row 0 in t) - T(row 1 in t)].  corner_mode 1 treats
// nodes on exactly one wall; 2 also treats nodes on two or more walls (corners, edges).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_wall_corr (int i, int j, int k, Array4<Real> const& rhs,
                             Array4<Real const> const& vel, Array4<int const> const& msk,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv, Box const& nddom,
                             GpuArray<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                             GpuArray<LinOpBCType,AMREX_SPACEDIM> const& hibc,
                             int corner_mode) noexcept
{
    if (msk(i,j,k)) { return; }
    IntVect const iv(AMREX_D_DECL(i,j,k));
    GpuArray<int,AMREX_SPACEDIM> wall;
    int const nwall = nodelap_detail::mlndlap_walls(iv, nddom, lobc, hibc, wall);
    if (nwall == 0 || (nwall > 1 && corner_mode < 2)) { return; }
    Real corr = Real(0.0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (int t = 0; t < AMREX_SPACEDIM; ++t) {
            if (t == d || wall[t] < 0) { continue; }
            corr += Real(0.25) * (nodelap_detail::mlndlap_wall_term(iv, d, t, 0, wall, vel, dxinv)
                                - nodelap_detail::mlndlap_wall_term(iv, d, t, 1, wall, vel, dxinv));
        }
    }
    rhs(i,j,k) += corr;
}

// The same for a cell-centred source: its control-volume mean is its value at h/4 from
// each wall, while mlndlap_rhcc (after the Neumann scaling) uses the first-row mean.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_rhcc_wall_corr (int i, int j, int k, Array4<Real> const& rhs,
                             Array4<Real const> const& rhcc, Array4<int const> const& msk,
                             Box const& nddom,
                             GpuArray<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                             GpuArray<LinOpBCType,AMREX_SPACEDIM> const& hibc,
                             int corner_mode) noexcept
{
    if (msk(i,j,k)) { return; }
    IntVect const iv(AMREX_D_DECL(i,j,k));
    GpuArray<int,AMREX_SPACEDIM> wall;
    int const nwall = nodelap_detail::mlndlap_walls(iv, nddom, lobc, hibc, wall);
    if (nwall == 0 || (nwall > 1 && corner_mode < 2)) { return; }
    Real corr = Real(0.0);
    for (int t = 0; t < AMREX_SPACEDIM; ++t) {
        if (wall[t] < 0) { continue; }
        corr += Real(0.25) * (nodelap_detail::mlndlap_wall_src(iv, t, 0, wall, rhcc)
                            - nodelap_detail::mlndlap_wall_src(iv, t, 1, wall, rhcc));
    }
    rhs(i,j,k) += corr;
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_sten (int i, int j, int k, Array4<T> const& x,
                             Array4<std::type_identity_t<T> const> const& sten,
                             Array4<int const> const& msk, std::type_identity_t<T> s0_norm0) noexcept
{
    if (!msk(i,j,k) && std::abs(sten(i,j,k,0)) > s0_norm0) {
        x(i,j,k) /= sten(i,j,k,0);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_sten (int i, int j, int k, Array4<T> const& sol,
                          std::type_identity_t<T> Ax, Array4<std::type_identity_t<T> const> const& rhs,
                          Array4<std::type_identity_t<T> const> const& sten,
                          Array4<int const> const& msk) noexcept
{
    if (msk(i,j,k)) {
        sol(i,j,k) = T(0.0);
    } else if (sten(i,j,k,0) != T(0.0)) {
        sol(i,j,k) += T(2./3.) * (rhs(i,j,k) - Ax) / sten(i,j,k,0);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_sten (Box const& bx, Array4<T> const& sol,
                          Array4<std::type_identity_t<T> const> const& Ax,
                          Array4<std::type_identity_t<T> const> const& rhs,
                          Array4<std::type_identity_t<T> const> const& sten,
                          Array4<int const> const& msk) noexcept
{
    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = T(0.0);
        } else if (sten(i,j,k,0) != T(0.0)) {
            sol(i,j,k) += T(2./3.) * (rhs(i,j,k) - Ax(i,j,k)) / sten(i,j,k,0);
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
        if (msk(i,j,k) == fine_flag) { return true; }
    }}}
    return false;
}

}

#endif
