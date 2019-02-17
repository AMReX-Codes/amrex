#ifndef AMREX_MLPOISSON_2D_K_H_
#define AMREX_MLPOISSON_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx (int i, int j, int /*k*/, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Real dhx, Real dhy,
                      Real const* AMREX_RESTRICT rc,
                      Real const* AMREX_RESTRICT re,
                      int rlo)
{
    y(i,j,0) = dhx * (re[i-rlo]*x(i-1,j,0)
                      - (re[i-rlo]+re[i+1-rlo])*x(i,j,0)
                      + re[i+1-rlo]*x(i+1,j,0))
        +   dhy * rc[i-rlo] * (x(i,j-1,0) - 2.*x(i,j,0) + x(i,j+1,0));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_x (Box const& box, Array4<Real> const& fx,
                       Array4<Real const> const& sol, Real dxinv,
                       Real const* AMREX_RESTRICT re, int rlo)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fx(i,j,0) = dxinv*re[i-rlo]*(sol(i,j,0)-sol(i-1,j,0));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_xface (Box const& box, Array4<Real> const& fx,
                           Array4<Real const> const& sol, Real dxinv, int xlen,
                           Real const* AMREX_RESTRICT re, int rlo)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        int i = lo.x;
        fx(i,j,0) = dxinv*re[i-rlo]*(sol(i,j,0)-sol(i-1,j,0));
        i += xlen;
        fx(i,j,0) = dxinv*re[i-rlo]*(sol(i,j,0)-sol(i-1,j,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_y (Box const& box, Array4<Real> const& fy,
                       Array4<Real const> const& sol, Real dyinv,
                       Real const* AMREX_RESTRICT rc, int rlo)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,0) = dyinv*rc[i-rlo]*(sol(i,j,0)-sol(i,j-1,0));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_yface (Box const& box, Array4<Real> const& fy,
                           Array4<Real const> const& sol, Real dyinv, int ylen,
                           Real const* AMREX_RESTRICT rc, int rlo)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int j = lo.y;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fy(i,j,0) = dyinv*rc[i-rlo]*(sol(i,j,0)-sol(i,j-1,0));
    }
    j += ylen;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fy(i,j,0) = dyinv*rc[i-rlo]*(sol(i,j,0)-sol(i,j-1,0));
    }
}

}

#endif
