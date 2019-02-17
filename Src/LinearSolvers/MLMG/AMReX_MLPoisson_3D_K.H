#ifndef AMREX_MLPOISSON_3D_K_H_
#define AMREX_MLPOISSON_3D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx (int i, int j, int k, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Real dhx, Real dhy, Real dhz)
{
    y(i,j,k) = dhx * (x(i-1,j,k) - 2.0*x(i,j,k) + x(i+1,j,k))
        +      dhy * (x(i,j-1,k) - 2.0*x(i,j,k) + x(i,j+1,k))
        +      dhz * (x(i,j,k-1) - 2.0*x(i,j,k) + x(i,j,k+1));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_x (Box const& box, Array4<Real> const& fx,
                       Array4<Real const> const& sol, Real dxinv)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fx(i,j,k) = dxinv*(sol(i,j,k)-sol(i-1,j,k));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_xface (Box const& box, Array4<Real> const& fx,
                           Array4<Real const> const& sol, Real dxinv, int xlen)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            int i = lo.x;
            fx(i,j,k) = dxinv*(sol(i,j,k)-sol(i-1,j,k));
            i += xlen;
            fx(i,j,k) = dxinv*(sol(i,j,k)-sol(i-1,j,k));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_y (Box const& box, Array4<Real> const& fy,
                       Array4<Real const> const& sol, Real dyinv)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fy(i,j,k) = dyinv*(sol(i,j,k)-sol(i,j-1,k));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_yface (Box const& box, Array4<Real> const& fy,
                           Array4<Real const> const& sol, Real dyinv, int ylen)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int k = lo.z; k <= hi.z; ++k) {
        int j = lo.y;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,k) = dyinv*(sol(i,j,k)-sol(i,j-1,k));
        }
        j += ylen;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,k) = dyinv*(sol(i,j,k)-sol(i,j-1,k));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_z (Box const& box, Array4<Real> const& fz,
                       Array4<Real const> const& sol, Real dzinv)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fz(i,j,k) = dzinv*(sol(i,j,k)-sol(i,j,k-1));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_zface (Box const& box, Array4<Real> const& fz,
                           Array4<Real const> const& sol, Real dzinv, int zlen)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = lo.z;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fz(i,j,k) = dzinv*(sol(i,j,k)-sol(i,j,k-1));
        }
    }

    k += zlen;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fz(i,j,k) = dzinv*(sol(i,j,k)-sol(i,j,k-1));
        }
    }
}

}

#endif
