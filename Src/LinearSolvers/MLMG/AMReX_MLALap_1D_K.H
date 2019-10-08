#ifndef AMREX_MLALAP_1D_K_H_
#define AMREX_MLALAP_1D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_adotx (Box const& box, Array4<Real> const& y,
                   Array4<Real const> const& x,
                   Array4<Real const> const& a,
                   GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                   Real alpha, Real beta) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        y(i,0,0) = alpha*a(i,0,0)*x(i,0,0)
            - dhx * (x(i+1,0,0) - 2.*x(i,0,0) + x(i-1,0,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_adotx_m (Box const& box, Array4<Real> const& y,
                     Array4<Real const> const& x,
                     Array4<Real const> const& a,
                     GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                     Real alpha, Real beta,
                     Real dx, Real probxlo) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        Real rel = (probxlo + i   *dx) * (probxlo + i   *dx);
        Real rer = (probxlo +(i+1)*dx) * (probxlo +(i+1)*dx);
        Real rc  = (probxlo +(i+.5)*dx) * (probxlo +(i+.5)*dx);
        y(i,0,0) = alpha*a(i,0,0)*x(i,0,0)*rc
            - dhx * (rer*(x(i+1,0,0) - x(i  ,0,0))
                   - rel*(x(i  ,0,0) - x(i-1,0,0)));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_normalize (Box const& box, Array4<Real> const& x,
                       Array4<Real const> const& a,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                       Real alpha, Real beta) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        x(i,0,0) /= alpha*a(i,0,0) + dhx*2.0;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_normalize_m (Box const& box, Array4<Real> const& x,
                         Array4<Real const> const& a,
                         GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                         Real alpha, Real beta, Real dx, Real probxlo) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        Real rel = (probxlo + i   *dx) * (probxlo + i   *dx);
        Real rer = (probxlo +(i+1)*dx) * (probxlo +(i+1)*dx);
        Real rc  = (probxlo +(i+.5)*dx) * (probxlo +(i+.5)*dx);
        x(i,0,0) /= alpha*a(i,0,0)*rc + dhx*(rel+rer);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_x (Box const& box, Array4<Real> const& fx,
                    Array4<Real const> const& sol, Real fac) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fx(i,0,0) = -fac*(sol(i,0,0)-sol(i-1,0,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_x_m (Box const& box, Array4<Real> const& fx,
                      Array4<Real const> const& sol, Real fac,
                      Real dx, Real probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        Real re = (probxlo + i*dx) * (probxlo + i*dx);
        fx(i,0,0) = -fac*re*(sol(i,0,0)-sol(i-1,0,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_xface (Box const& box, Array4<Real> const& fx,
                        Array4<Real const> const& sol, Real fac, int xlen) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int i = lo.x;
    fx(i,0,0) = -fac*(sol(i,0,0)-sol(i-1,0,0));
    i += xlen;
    fx(i,0,0) = -fac*(sol(i,0,0)-sol(i-1,0,0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_xface_m (Box const& box, Array4<Real> const& fx,
                          Array4<Real const> const& sol, Real fac, int xlen,
                          Real dx, Real probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int i = lo.x;
    Real re = (probxlo + i*dx) * (probxlo + i*dx);
    fx(i,0,0) = -fac*re*(sol(i,0,0)-sol(i-1,0,0));
    i += xlen;
    re = (probxlo + i*dx) * (probxlo + i*dx);
    fx(i,0,0) = -fac*re*(sol(i,0,0)-sol(i-1,0,0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_gsrb (Box const& box, Array4<Real> const& phi,
                  Array4<Real const> const& rhs, Real alpha,
                  Real dhx, Array4<Real const> const& a,
                  Array4<Real const> const& f0, Array4<int const> const& m0,
                  Array4<Real const> const& f1, Array4<int const> const& m1,
                  Box const& vbox, int redblack) noexcept
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
            
            Real delta = dhx*(cf0+cf1);
            
            Real gamma = alpha*a(i,0,0) + dhx*2.0;
            
            Real rho = dhx*(phi(i-1,0,0) + phi(i+1,0,0));

            phi(i,0,0) = (rhs(i,0,0) + rho - phi(i,0,0)*delta)
                / (gamma - delta);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_gsrb_m (Box const& box, Array4<Real> const& phi,
                    Array4<Real const> const& rhs, Real alpha,
                    Real dhx, Array4<Real const> const& a,
                    Array4<Real const> const& f0, Array4<int const> const& m0,
                    Array4<Real const> const& f1, Array4<int const> const& m1,
                    Box const& vbox, int redblack,
                    Real dx, Real probxlo)
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

            Real rel = (probxlo + i   *dx) * (probxlo + i   *dx);
            Real rer = (probxlo +(i+1)*dx) * (probxlo +(i+1)*dx);
            Real rc  = (probxlo +(i+.5)*dx) * (probxlo +(i+.5)*dx);
            
            Real delta = dhx*(rel*cf0 + rer*cf1);
            
            Real gamma = alpha*a(i,0,0)*rc + dhx*(rel+rer);
            
            Real rho = dhx*(rel*phi(i-1,0,0) + rer*phi(i+1,0,0));

            phi(i,0,0) = (rhs(i,0,0) + rho - phi(i,0,0)*delta)
                / (gamma - delta);
        }
    }
}

}
#endif
