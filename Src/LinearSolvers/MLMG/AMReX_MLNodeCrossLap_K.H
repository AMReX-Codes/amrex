#ifndef AMREX_MLNODECROSSLAP_K_H_
#define AMREX_MLNODECROSSLAP_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_Box.H>

#include <limits>

#if (AMREX_SPACEDIM == 2)

namespace amrex {

// Stencil of D sigma G: corners cc, x neighbors cx, y neighbors -cx, center c0.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndcrosslap_adotx (int i, int j, int k, Array4<Real const> const& x,
                         Real cc, Real cx, Real c0) noexcept
{
    return cc * (x(i-1,j-1,k) + x(i+1,j-1,k) + x(i-1,j+1,k) + x(i+1,j+1,k))
        +  cx * (x(i-1,j,k) + x(i+1,j,k) - x(i,j-1,k) - x(i,j+1,k))
        +  c0 * x(i,j,k);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndcrosslap_gauss_seidel (int i, int j, int k, Array4<Real> const& sol,
                                Array4<Real const> const& rhs, Array4<int const> const& msk,
                                Real cc, Real cx, Real c0) noexcept
{
    if (msk(i,j,k)) {
        sol(i,j,k) = Real(0.0);
    } else {
        Real const Ax = mlndcrosslap_adotx(i, j, k, sol, cc, cx, c0);
        sol(i,j,k) += (rhs(i,j,k) - Ax) / c0;
    }
}

// Matrix row of node (i,j) without the x and y neighbor couplings, from the
// cells inside ccdom (grown at periodic boundaries). Slot 0 is the diagonal
// and slots 1-4 are the corners. Returns whether slot `slot` has an entry,
// and sets col/val for it.
template <typename AlgGid>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool mlndcrosslap_matrix_entry (int i, int j, int k, int slot,
                                Array4<AlgGid const> const& gid, Box const& ccdom,
                                Real cc, AlgGid& col, Real& val) noexcept
{
    constexpr auto gidmax = std::numeric_limits<AlgGid>::max();
    // Corner node and the cell between it and (i,j).
    constexpr int di[4] = {-1, 1, -1, 1};
    constexpr int dj[4] = {-1, -1, 1, 1};
    if (slot == 0) {
        Real diag = Real(0.0);
        for (int m = 0; m < 4; ++m) {
            if (ccdom.contains(i + (di[m]-1)/2, j + (dj[m]-1)/2, k)) { diag -= cc; }
        }
        col = gid(i,j,k);
        val = diag;
        return true;
    } else {
        int const m = slot - 1;
        if (ccdom.contains(i + (di[m]-1)/2, j + (dj[m]-1)/2, k) &&
            gid(i+di[m],j+dj[m],k) < gidmax) {
            col = gid(i+di[m],j+dj[m],k);
            val = cc;
            return true;
        } else {
            return false;
        }
    }
}

}

#endif
#endif
