#ifndef AMREX_MLABECLAP_1D_K_H_
#define AMREX_MLABECLAP_1D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_adotx (Box const& box, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Array4<Real const> const& a,
                      Array4<Real const> const& bX,
                      GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                      Real alpha, Real beta)
{
    const Real dhx = beta*dxinv[0]*dxinv[0];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        y(i,0,0) = alpha*a(i,0,0)*x(i,0,0)
            - dhx * (bX(i+1,0,0)*(x(i+1,0,0) - x(i  ,0,0))
                   - bX(i  ,0,0)*(x(i  ,0,0) - x(i-1,0,0)));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_normalize (Box const& box, Array4<Real> const& x,
                          Array4<Real const> const& a,
                          Array4<Real const> const& bX,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                          Real alpha, Real beta)
{
    const Real dhx = beta*dxinv[0]*dxinv[0];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        x(i,0,0) /= alpha*a(i,0,0) + dhx*(bX(i,0,0)+bX(i+1,0,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_flux_x (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                       Array4<Real const> const& bx, Real fac)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fx(i,0,0) = -fac*bx(i,0,0)*(sol(i,0,0)-sol(i-1,0,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_flux_xface (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                           Array4<Real const> const& bx, Real fac, int xlen)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int i = lo.x;
    fx(i,0,0) = -fac*bx(i,0,0)*(sol(i,0,0)-sol(i-1,0,0));
    i += xlen;
    fx(i,0,0) = -fac*bx(i,0,0)*(sol(i,0,0)-sol(i-1,0,0));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void abec_gsrb (Box const& box, Array4<Real> const& phi,
                Array4<Real const> const& rhs, Real alpha,
                Real dhx, Array4<Real const> const& a,
                Array4<Real const> const& bX,
                Array4<Real const> const& f0, Array4<int const> const& m0,
                Array4<Real const> const& f1, Array4<int const> const& m1,
                Box const& vbox, int nc, int redblack)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    for (int n = 0; n < nc; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((i+redblack)%2 == 0) {
                Real cf0 = (i == vlo.x and m0(vlo.x-1,0,0) > 0)
                    ? f0(vlo.x,0,0) : 0.0;
                Real cf1 = (i == vhi.x and m1(vhi.x+1,0,0) > 0)
                    ? f1(vhi.x,0,0) : 0.0;

                Real delta = dhx*(bX(i,0,0)*cf0 + bX(i+1,0,0)*cf1);

                Real gamma = alpha*a(i,0,0)
                    +   dhx*( bX(i,0,0) + bX(i+1,0,0) );

                Real rho = dhx*(bX(i  ,0  ,0)*phi(i-1,0  ,0,n)
                              + bX(i+1,0  ,0)*phi(i+1,0  ,0,n));

                phi(i,0,0,n) = (rhs(i,0,0,n) + rho - phi(i,0,0,n)*delta)
                    / (gamma - delta);
            }
        }
    }
}

}
#endif
