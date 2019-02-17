#ifndef AMREX_MLPOISSON_2D_K_H_
#define AMREX_MLPOISSON_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx (int i, int /*j*/, int /*k*/, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Real dhx,
                      Real const* AMREX_RESTRICT re,
                      int rlo)
{
    y(i,0,0) = dhx * (re[i-rlo]*x(i-1,0,0)
                      - (re[i-rlo]+re[i+1-rlo])*x(i,0,0)
                      + re[i+1-rlo]*x(i+1,0,0));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_x (Box const& box, Array4<Real> const& fx,
                       Array4<Real const> const& sol, Real dxinv,
                       Real const* AMREX_RESTRICT re, int rlo)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fx(i,0,0) = dxinv*re[i-rlo]*(sol(i,0,0)-sol(i-1,0,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_xface (Box const& box, Array4<Real> const& fx,
                           Array4<Real const> const& sol, Real dxinv, int xlen,
                           Real const* AMREX_RESTRICT re, int rlo)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int i = lo.x;
    fx(i,0,0) = dxinv*re[i-rlo]*(sol(i,0,0)-sol(i-1,0,0));
    i += xlen;
    fx(i,0,0) = dxinv*re[i-rlo]*(sol(i,0,0)-sol(i-1,0,0));
}

}

#endif
