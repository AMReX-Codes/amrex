#ifndef AMREX_MLPOISSON_2D_K_H_
#define AMREX_MLPOISSON_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx (int i, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Real dhx) noexcept
{
    y(i,0,0) = dhx * (x(i-1,0,0) - 2.0*x(i,0,0) + x(i+1,0,0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx_m (int i, Array4<Real> const& y,
                        Array4<Real const> const& x,
                        Real dhx, Real dx, Real probxlo) noexcept
{
    Real rel = (probxlo + i   *dx) * (probxlo + i   *dx);
    Real rer = (probxlo +(i+1)*dx) * (probxlo +(i+1)*dx);
    y(i,0,0) = dhx * (rel*x(i-1,0,0) - (rel+rer)*x(i,0,0) + rer*x(i+1,0,0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_x (Box const& box, Array4<Real> const& fx,
                       Array4<Real const> const& sol, Real dxinv) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fx(i,0,0) = dxinv*(sol(i,0,0)-sol(i-1,0,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_x_m (Box const& box, Array4<Real> const& fx,
                         Array4<Real const> const& sol, Real dxinv,
                         Real dx, Real probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        Real re = (probxlo + i*dx) * (probxlo + i*dx);
        fx(i,0,0) = dxinv*re*(sol(i,0,0)-sol(i-1,0,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_xface (Box const& box, Array4<Real> const& fx,
                           Array4<Real const> const& sol, Real dxinv, int xlen) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int i = lo.x;
    fx(i,0,0) = dxinv*(sol(i,0,0)-sol(i-1,0,0));
    i += xlen;
    fx(i,0,0) = dxinv*(sol(i,0,0)-sol(i-1,0,0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_xface_m (Box const& box, Array4<Real> const& fx,
                             Array4<Real const> const& sol, Real dxinv, int xlen,
                             Real dx, Real probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int i = lo.x;
    Real re = (probxlo + i*dx) * (probxlo + i*dx);    
    fx(i,0,0) = dxinv*re*(sol(i,0,0)-sol(i-1,0,0));
    i += xlen;
    re = (probxlo + i*dx) * (probxlo + i*dx);
    fx(i,0,0) = dxinv*re*(sol(i,0,0)-sol(i-1,0,0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_gsrb (Box const& box, Array4<Real> const& phi, Array4<Real const> const& rhs,
                     Real dhx,
                     Array4<Real const> const& f0, Array4<int const> const& m0,
                     Array4<Real const> const& f1, Array4<int const> const& m1,
                     Box const& vbox, int redblack) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    Real gamma = -dhx*2.0;

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        if ((i+redblack)%2 == 0) {
            Real cf0 = (i == vlo.x and m0(vlo.x-1,0,0) > 0)
                ? f0(vlo.x,0,0) : 0.0;
            Real cf1 = (i == vhi.x and m1(vhi.x+1,0,0) > 0)
                ? f1(vhi.x,0,0) : 0.0;

            Real g_m_d = gamma + dhx*(cf0+cf1);

            Real res = rhs(i,0,0) - gamma*phi(i,0,0)
                - dhx*(phi(i-1,0,0) + phi(i+1,0,0));

            phi(i,0,0) = phi(i,0,0) + res /g_m_d;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_gsrb_m (Box const& box, Array4<Real> const& phi, Array4<Real const> const& rhs,
                       Real dhx,
                       Array4<Real const> const& f0, Array4<int const> const& m0,
                       Array4<Real const> const& f1, Array4<int const> const& m1,
                       Box const& vbox, int redblack, Real dx, Real probxlo) noexcept
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

            Real gamma = -dhx*(rel+rer);

            Real g_m_d = gamma + dhx*(rel*cf0+rer*cf1);

            Real res = rhs(i,0,0) - gamma*phi(i,0,0)
                - dhx*(rel*phi(i-1,0,0) + rer*phi(i+1,0,0));
       
            phi(i,0,0) = phi(i,0,0) + res /g_m_d;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_normalize (Box const& box, Array4<Real> const& x,
                          Real dhx, Real dx, Real probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        Real rel = (probxlo + i   *dx) * (probxlo + i   *dx);
        Real rer = (probxlo +(i+1)*dx) * (probxlo +(i+1)*dx);
        x(i,0,0) /= (-dhx*(rel+rer));
    }
}

}

#endif
