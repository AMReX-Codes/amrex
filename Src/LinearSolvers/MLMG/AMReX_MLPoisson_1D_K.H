#ifndef AMREX_MLPOISSON_2D_K_H_
#define AMREX_MLPOISSON_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx (int i, int /*j*/, int /*k*/, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Real dhx,
                      Real const* AMREX_RESTRICT re,
                      int rlo) noexcept
{
    y(i,0,0) = dhx * (re[i-rlo]*x(i-1,0,0)
                      - (re[i-rlo]+re[i+1-rlo])*x(i,0,0)
                      + re[i+1-rlo]*x(i+1,0,0));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_x (Box const& box, Array4<Real> const& fx,
                       Array4<Real const> const& sol, Real dxinv,
                       Real const* AMREX_RESTRICT re, int rlo) noexcept
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
                           Real const* AMREX_RESTRICT re, int rlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int i = lo.x;
    fx(i,0,0) = dxinv*re[i-rlo]*(sol(i,0,0)-sol(i-1,0,0));
    i += xlen;
    fx(i,0,0) = dxinv*re[i-rlo]*(sol(i,0,0)-sol(i-1,0,0));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_gsrb (Box const& box, Array4<Real> const& phi, Array4<Real const> const& rhs,
                     Real dhx,
                     Array4<Real const> const& f0, Array4<int const> const& m0,
                     Array4<Real const> const& f1, Array4<int const> const& m1,
                     Box const& vbox, int redblack,
                     Real const* AMREX_RESTRICT re, int rlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        if ((i+redblack)%2 == 0) {
            Real cf0 = (i == vlo.x and m0(vlo.x-1,0,0) > 0)
                ? f0(vlo.x,0,0) : 0.0;
            Real cf1 = (i == vhi.x and m1(vhi.x+1,0,0) > 0)
                ? f1(vhi.x,0,0) : 0.0;

            Real gamma = -dhx*(re[i-rlo]+re[i+1-rlo]);

            Real g_m_d = gamma + dhx*(re[i-rlo]*cf0+re[i+1-rlo]*cf1);

            Real res = rhs(i,0,0) - gamma*phi(i,0,0)
                - dhx*(re[i-rlo]*phi(i-1,0,0) + re[i+1-rlo]*phi(i+1,0,0));
       
            phi(i,0,0) = phi(i,0,0) + res /g_m_d;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_normalize (Box const& box, Array4<Real> const& x,
                          Real const* AMREX_RESTRICT re, int rlo, Real dhx) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        x(i,0,0) /= (-dhx*(re[i-rlo]+re[i+1-rlo]));
    }
}

}

#endif
