#ifndef AMREX_MLPOISSON_2D_K_H_
#define AMREX_MLPOISSON_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx (int i, int j, int /*k*/, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Real dhx, Real dhy,
                      Real const* AMREX_RESTRICT rc,
                      Real const* AMREX_RESTRICT re,
                      int rlo) noexcept
{
    y(i,j,0) = dhx * (re[i-rlo]*x(i-1,j,0)
                      - (re[i-rlo]+re[i+1-rlo])*x(i,j,0)
                      + re[i+1-rlo]*x(i+1,j,0))
        +   dhy * rc[i-rlo] * (x(i,j-1,0) - 2.*x(i,j,0) + x(i,j+1,0));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_flux_x (Box const& box, Array4<Real> const& fx,
                       Array4<Real const> const& sol, Real dxinv,
                       Real const* AMREX_RESTRICT re, int rlo) noexcept
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
                           Real const* AMREX_RESTRICT re, int rlo) noexcept
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
                       Real const* AMREX_RESTRICT rc, int rlo) noexcept
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
                           Real const* AMREX_RESTRICT rc, int rlo) noexcept
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

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_gsrb (Box const& box, Array4<Real> const& phi, Array4<Real const> const& rhs,
                     Real dhx, Real dhy,
                     Array4<Real const> const& f0, Array4<int const> const& m0,
                     Array4<Real const> const& f1, Array4<int const> const& m1,
                     Array4<Real const> const& f2, Array4<int const> const& m2,
                     Array4<Real const> const& f3, Array4<int const> const& m3,
                     Box const& vbox, int redblack,
                     Real const* AMREX_RESTRICT rc, Real const* AMREX_RESTRICT re, int rlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((i+j+redblack)%2 == 0) {
                Real cf0 = (i == vlo.x and m0(vlo.x-1,j,0) > 0)
                    ? f0(vlo.x,j,0) : 0.0;
                Real cf1 = (j == vlo.y and m1(i,vlo.y-1,0) > 0)
                    ? f1(i,vlo.y,0) : 0.0;
                Real cf2 = (i == vhi.x and m2(vhi.x+1,j,0) > 0)
                    ? f2(vhi.x,j,0) : 0.0;
                Real cf3 = (j == vhi.y and m3(i,vhi.y+1,0) > 0)
                    ? f3(i,vhi.y,0) : 0.0;

                Real gamma = -dhx*(re[i-rlo]+re[i+1-rlo]) - 2.0*dhy*rc[i-rlo];

                Real g_m_d = gamma + dhx*(re[i-rlo]*cf0+re[i+1-rlo]*cf2) + dhy*rc[i-rlo]*(cf1+cf3);

                Real res = rhs(i,j,0) - gamma*phi(i,j,0)
                    - dhx*(re[i-rlo]*phi(i-1,j,0) + re[i+1-rlo]*phi(i+1,j,0))
                    - dhy*rc[i-rlo]*(phi(i,j-1,0) + phi(i,j+1,0));

                phi(i,j,0) = phi(i,j,0) + res /g_m_d;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_normalize (Box const& box, Array4<Real> const& x,
                          Real const* AMREX_RESTRICT rc, Real const* AMREX_RESTRICT re, int rlo,
                          Real dhx, Real dhy) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            x(i,j,0) /= (-dhx*(re[i-rlo]+re[i+1-rlo]) - dhy*rc[i-rlo]*2.0);
        }
    }
}

}

#endif
