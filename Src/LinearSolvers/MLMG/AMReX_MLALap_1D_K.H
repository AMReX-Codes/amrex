#ifndef AMREX_MLALAP_1D_K_H_
#define AMREX_MLALAP_1D_K_H_

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

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        y(i,0,0) = alpha*a(i,0,0)*x(i,0,0)*rc[i-rlo]
            - dhx * (re[i+1-rlo]*(x(i+1,0,0) - x(i  ,0,0))
                   - re[i  -rlo]*(x(i  ,0,0) - x(i-1,0,0)));
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

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        x(i,0,0) /= alpha*a(i,0,0)*rc[i-rlo] + dhx*(re[i-rlo]+re[i+1-rlo]);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_flux_x (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                    Real fac, Real const* AMREX_RESTRICT re, int rlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fx(i,0,0) = -fac*re[i-rlo]*(sol(i,0,0)-sol(i-1,0,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_flux_xface (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                        Real fac, int xlen, Real const* AMREX_RESTRICT re, int rlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int i = lo.x;
    fx(i,0,0) = -fac*re[i-rlo]*(sol(i,0,0)-sol(i-1,0,0));
    i += xlen;
    fx(i,0,0) = -fac*re[i-rlo]*(sol(i,0,0)-sol(i-1,0,0));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_gsrb (Box const& box, Array4<Real> const& phi,
                  Array4<Real const> const& rhs, Real alpha,
                  Real dhx, Array4<Real const> const& a,
                  Array4<Real const> const& f0, Array4<int const> const& m0,
                  Array4<Real const> const& f1, Array4<int const> const& m1,
                  Box const& vbox, int redblack,
                  Real const* AMREX_RESTRICT rc,
                  Real const* AMREX_RESTRICT re,
                  int rlo) noexcept
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
            
            Real delta = dhx*(re[i-rlo]*cf0 + re[i+1-rlo]*cf1);
            
            Real gamma = alpha*a(i,0,0)*rc[i-rlo]
                +   dhx*(re[i-rlo]+re[i+1-rlo]);
            
            Real rho = dhx*(re[i  -rlo]*phi(i-1,0,0)
                          + re[i+1-rlo]*phi(i+1,0,0));

            phi(i,0,0) = (rhs(i,0,0) + rho - phi(i,0,0)*delta)
                / (gamma - delta);
        }
    }
}

}
#endif
