#ifndef AMREX_MLABECLAP_DRIFT_K_H_
#define AMREX_MLABECLAP_DRIFT_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>

/**
 * \file AMReX_MLABecLap_Drift_K.H
 *
 * Kernels of MLABecLaplacian with the upwind drift term of setDriftCoeffs,
 *
 *   alpha a phi - beta div(b grad phi) - cscalar div(c phi),
 *
 * where c are face velocities and the flux through a face is c times the cell
 * value upwind of it. For one direction, the row of cell i is
 *
 *   diag phi_i - wm phi_{i-1} - wp phi_{i+1},
 *
 * with the diffusive and drift parts of the face weights wm and wp.
 */

namespace amrex {

//! The coefficients of the operator with drift on one box
template <typename T>
struct MLABecLapDriftCoeffs
{
    Array4<T const> a;
    GpuArray<Array4<T const>,AMREX_SPACEDIM> b, c;
    GpuArray<T,AMREX_SPACEDIM> dxinv;
    T alpha, beta, cscalar;
};

//! The coefficients of the operator with drift on all the boxes of a level
template <typename T>
struct MLABecLapDriftCoeffsMA
{
    MultiArray4<T const> a;
    GpuArray<MultiArray4<T const>,AMREX_SPACEDIM> b, c;
    GpuArray<T,AMREX_SPACEDIM> dxinv;
    T alpha, beta, cscalar;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    MLABecLapDriftCoeffs<T> operator[] (int box_no) const noexcept
    {
        MLABecLapDriftCoeffs<T> r{a[box_no], {}, {}, dxinv, alpha, beta, cscalar};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            r.b[d] = b[d][box_no];
            r.c[d] = c[d][box_no];
        }
        return r;
    }
};

/**
 * \brief The row of cell (iv,n): returns the diagonal, and sets \p off to the
 * sum of the neighbor weights times \p x (when \p x is given) and \p delta to
 * the boundary folding of mlabeclap_relax_drift (when \p m is).
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T mlabeclap_drift_row (IntVect const& iv, int n, MLABecLapDriftCoeffs<T> const& op,
                       Array4<T const> const& x, T& off, T& delta,
                       GpuArray<Array4<int const>,2*AMREX_SPACEDIM> const& m = {},
                       GpuArray<Array4<T const>,2*AMREX_SPACEDIM> const& f = {},
                       Box const& vbox = Box()) noexcept
{
    T diag = op.alpha*op.a(iv);
    off = T(0.0);
    delta = T(0.0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        IntVect const e = IntVect::TheDimensionVector(d);
        T const dhb = op.beta*op.dxinv[d]*op.dxinv[d];
        T const dhc = op.cscalar*op.dxinv[d];
        T const bm = op.b[d](iv,n);
        T const bp = op.b[d](iv+e,n);
        T const um = op.c[d](iv,n);
        T const up = op.c[d](iv+e,n);
        diag += dhb*(bm + bp);
        T wm = dhb*bm;
        T wp = dhb*bp;
        // The low face adds +dhc*um*phi_upwind, the high face -dhc*up*phi_upwind
        if (um > T(0.0)) { wm -= dhc*um; } else { diag += dhc*um; }
        if (up > T(0.0)) { diag -= dhc*up; } else { wp += dhc*up; }
        if (x) { off += wm*x(iv-e,n) + wp*x(iv+e,n); }
        if (m[d]) {
            if (iv[d] == vbox.smallEnd(d) && m[d](iv-e) > 0) {
                delta += wm*f[d](iv,n);
            }
            if (iv[d] == vbox.bigEnd(d) && m[d+AMREX_SPACEDIM](iv+e) > 0) {
                delta += wp*f[d+AMREX_SPACEDIM](iv,n);
            }
        }
    }
    return diag;
}

//! y = A x at (i,j,k,n)
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_adotx_drift (int i, int j, [[maybe_unused]] int k, int n, Array4<T> const& y,
                            Array4<T const> const& x, MLABecLapDriftCoeffs<T> const& op) noexcept
{
    IntVect const iv(AMREX_D_DECL(i,j,k));
    T off, delta;
    T const diag = mlabeclap_drift_row(iv, n, op, x, off, delta);
    y(iv,n) = diag*x(iv,n) - off;
}

//! x /= the diagonal of A at (i,j,k,n)
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_normalize_drift (int i, int j, [[maybe_unused]] int k, int n, Array4<T> const& x,
                                MLABecLapDriftCoeffs<T> const& op) noexcept
{
    IntVect const iv(AMREX_D_DECL(i,j,k));
    T off, delta;
    x(iv,n) /= mlabeclap_drift_row(iv, n, op, Array4<T const>{}, off, delta);
}

/**
 * \brief Relaxation of cell (i,j,k,n): phi += omega * (rhs - A phi) / D.
 *
 * D is the diagonal of A with the boundary conditions: where the ghost cell
 * across a face of the valid box \p vbox is set from the cells inside (mask
 * \p m > 0), its weight times the coefficient \p f of this cell in it is taken
 * off the diagonal. With \p Ax, a Jacobi update with A phi given; without, a
 * Gauss-Seidel update (the caller visits the cells of one color only).
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_relax_drift (int i, int j, [[maybe_unused]] int k, int n, Array4<T> const& phi,
                            Array4<T const> const& rhs, Array4<T const> const& Ax,
                            MLABecLapDriftCoeffs<T> const& op,
                            GpuArray<Array4<int const>,2*AMREX_SPACEDIM> const& m,
                            GpuArray<Array4<T const>,2*AMREX_SPACEDIM> const& f,
                            T omega, Box const& vbox) noexcept
{
    IntVect const iv(AMREX_D_DECL(i,j,k));
    T off, delta;
    T const diag = mlabeclap_drift_row(iv, n, op, Ax ? Array4<T const>{} : Array4<T const>(phi),
                                       off, delta, m, f, vbox);
    T const res = Ax ? rhs(iv,n) - Ax(iv,n) : rhs(iv,n) - (diag*phi(iv,n) - off);
    // The sign of a stable diagonal; inflow through a Neumann face, or a cell
    // that only takes in with A = B = 0, breaks it
    AMREX_ASSERT((diag - delta) * ((op.beta != T(0.0)) ? op.beta : -op.cscalar) > T(0.0));
    phi(iv,n) += omega * res / (diag - delta);
}

/**
 * \brief Add the drift flux -cscalar * c * phi_upwind to \p flux on face
 * \p iv of direction \p d.
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_drift (IntVect const& iv, int d, int n, Array4<T> const& flux,
                           Array4<T const> const& sol, Array4<T const> const& c,
                           T cscalar) noexcept
{
    T const u = c(iv,n);
    T const s = (u > T(0.0)) ? sol(iv-IntVect::TheDimensionVector(d),n) : sol(iv,n);
    flux(iv,n) -= cscalar*u*s;
}

}

#endif
