#ifndef AMREX_MLALAP_2D_K_H_
#define AMREX_MLALAP_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_adotx (Box const& box, Array4<Real> const& y,
                   Array4<Real const> const& x,
                   Array4<Real const> const& a,
                   GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                   Real alpha, Real beta,
                   Real const* AMREX_RESTRICT rc,
                   Real const* AMREX_RESTRICT re,
                   int rlo) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            y(i,j,0) = alpha*a(i,j,0)*x(i,j,0)*rc[i-rlo]
                - dhx * (re[i+1-rlo]*(x(i+1,j,0) - x(i  ,j,0))
                       - re[i  -rlo]*(x(i  ,j,0) - x(i-1,j,0)))
                - dhy * rc[i-rlo] * (x(i,j-1,0) - 2.*x(i,j,0) + x(i,j+1,0));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_normalize (Box const& box, Array4<Real> const& x,
                       Array4<Real const> const& a,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                       Real alpha, Real beta,
                       Real const* AMREX_RESTRICT rc,
                       Real const* AMREX_RESTRICT re,
                       int rlo) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            x(i,j,0) /= alpha*a(i,j,0)*rc[i-rlo]
                + dhx*(re[i-rlo]+re[i+1-rlo])
                + dhy*(rc[i-rlo]*2.0);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_flux_x (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                    Real fac, Real const* AMREX_RESTRICT re, int rlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fx(i,j,0) = -fac*re[i-rlo]*(sol(i,j,0)-sol(i-1,j,0));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_flux_xface (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                        Real fac, int xlen, Real const* AMREX_RESTRICT re, int rlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        int i = lo.x;
        fx(i,j,0) = -fac*re[i-rlo]*(sol(i,j,0)-sol(i-1,j,0));
        i += xlen;
        fx(i,j,0) = -fac*re[i-rlo]*(sol(i,j,0)-sol(i-1,j,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_flux_y (Box const& box, Array4<Real> const& fy, Array4<Real const> const& sol,
                    Real fac, Real const* AMREX_RESTRICT rc, int rlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,0) = -fac*rc[i-rlo]*(sol(i,j,0)-sol(i,j-1,0));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_flux_yface (Box const& box, Array4<Real> const& fy, Array4<Real const> const& sol,
                        Real fac, int ylen, Real const* AMREX_RESTRICT rc, int rlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int j = lo.y;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fy(i,j,0) = -fac*rc[i-rlo]*(sol(i,j,0)-sol(i,j-1,0));
    }
    j += ylen;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fy(i,j,0) = -fac*rc[i-rlo]*(sol(i,j,0)-sol(i,j-1,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_gsrb (Box const& box, Array4<Real> const& phi,
                  Array4<Real const> const& rhs, Real alpha,
                  Real dhx, Real dhy, Array4<Real const> const& a,
                  Array4<Real const> const& f0, Array4<int const> const& m0,
                  Array4<Real const> const& f1, Array4<int const> const& m1,
                  Array4<Real const> const& f2, Array4<int const> const& m2,
                  Array4<Real const> const& f3, Array4<int const> const& m3,
                  Box const& vbox, int redblack,
                  Real const* AMREX_RESTRICT rc,
                  Real const* AMREX_RESTRICT re,
                  int rlo) noexcept
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

                Real delta = dhx*(re[i-rlo]*cf0 + re[i+1-rlo]*cf2)
                          +  dhy*rc[i-rlo]*(cf1 + cf3);

                Real gamma = alpha*a(i,j,0)*rc[i-rlo]
                    +   dhx*( re[i-rlo] + re[i+1-rlo] )
                    +   dhy*(2.*rc[i-rlo]);

                Real rho = dhx*(re[i  -rlo]*phi(i-1,j,0)
                              + re[i+1-rlo]*phi(i+1,j,0))
                    +dhy*rc[i-rlo]*(phi(i,j-1,0) + phi(i,j+1,0));

                phi(i,j,0) = (rhs(i,j,0) + rho - phi(i,j,0)*delta)
                    / (gamma - delta);
            }
        }
    }
}

}
#endif
