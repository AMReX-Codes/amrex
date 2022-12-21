#ifndef AMREX_MLPOISSON_2D_K_H_
#define AMREX_MLPOISSON_2D_K_H_
#include <AMReX_Config.H>

namespace amrex {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx (int i, Array4<T> const& y,
                      Array4<T const> const& x,
                      T dhx) noexcept
{
    y(i,0,0) = dhx * (x(i-1,0,0) - T(2.0)*x(i,0,0) + x(i+1,0,0));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx_os (int i, Array4<T> const& y,
                         Array4<T const> const& x,
                         Array4<int const> const& osm,
                         T dhx) noexcept
{
    if (osm(i,0,0) == 0) {
        y(i,0,0) = T(0.0);
    } else {
        y(i,0,0) = dhx * (x(i-1,0,0) - T(2.0)*x(i,0,0) + x(i+1,0,0));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx_m (int i, Array4<T> const& y,
                        Array4<T const> const& x,
                        T dhx, T dx, T probxlo) noexcept
{
    T rel = (probxlo + i   *dx) * (probxlo + i   *dx);
    T rer = (probxlo +(i+1)*dx) * (probxlo +(i+1)*dx);
    y(i,0,0) = dhx * (rel*x(i-1,0,0) - (rel+rer)*x(i,0,0) + rer*x(i+1,0,0));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_x (Box const& box, Array4<T> const& fx,
                       Array4<T const> const& sol, T dxinv) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fx(i,0,0) = dxinv*(sol(i,0,0)-sol(i-1,0,0));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_x_m (Box const& box, Array4<T> const& fx,
                         Array4<T const> const& sol, T dxinv,
                         T dx, T probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        T re = (probxlo + i*dx) * (probxlo + i*dx);
        fx(i,0,0) = dxinv*re*(sol(i,0,0)-sol(i-1,0,0));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_xface (Box const& box, Array4<T> const& fx,
                           Array4<T const> const& sol, T dxinv, int xlen) noexcept
{
    const auto lo = amrex::lbound(box);

    int i = lo.x;
    fx(i,0,0) = dxinv*(sol(i,0,0)-sol(i-1,0,0));
    i += xlen;
    fx(i,0,0) = dxinv*(sol(i,0,0)-sol(i-1,0,0));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_xface_m (Box const& box, Array4<T> const& fx,
                             Array4<T const> const& sol, T dxinv, int xlen,
                             T dx, T probxlo) noexcept
{
    const auto lo = amrex::lbound(box);

    int i = lo.x;
    T re = (probxlo + i*dx) * (probxlo + i*dx);
    fx(i,0,0) = dxinv*re*(sol(i,0,0)-sol(i-1,0,0));
    i += xlen;
    re = (probxlo + i*dx) * (probxlo + i*dx);
    fx(i,0,0) = dxinv*re*(sol(i,0,0)-sol(i-1,0,0));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_gsrb (int i, int, int, Array4<T> const& phi, Array4<T const> const& rhs,
                     T dhx,
                     Array4<T const> const& f0, Array4<int const> const& m0,
                     Array4<T const> const& f1, Array4<int const> const& m1,
                     Box const& vbox, int redblack) noexcept
{
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    T gamma = -dhx*T(2.0);

    if ((i+redblack)%2 == 0) {
        T cf0 = (i == vlo.x && m0(vlo.x-1,0,0) > 0)
            ? f0(vlo.x,0,0) : T(0.0);
        T cf1 = (i == vhi.x && m1(vhi.x+1,0,0) > 0)
            ? f1(vhi.x,0,0) : T(0.0);

        T g_m_d = gamma + dhx*(cf0+cf1);

        T res = rhs(i,0,0) - gamma*phi(i,0,0)
            - dhx*(phi(i-1,0,0) + phi(i+1,0,0));

        phi(i,0,0) = phi(i,0,0) + res /g_m_d;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_gsrb_os (int i, int, int, Array4<T> const& phi, Array4<T const> const& rhs,
                        Array4<int const> const& osm, T dhx,
                        Array4<T const> const& f0, Array4<int const> const& m0,
                        Array4<T const> const& f1, Array4<int const> const& m1,
                        Box const& vbox, int redblack) noexcept
{
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    T gamma = -dhx*T(2.0);

    if ((i+redblack)%2 == 0) {
        if (osm(i,0,0) == 0) {
            phi(i,0,0) = T(0.0);
        } else {
            T cf0 = (i == vlo.x && m0(vlo.x-1,0,0) > 0)
                ? f0(vlo.x,0,0) : T(0.0);
            T cf1 = (i == vhi.x && m1(vhi.x+1,0,0) > 0)
                ? f1(vhi.x,0,0) : T(0.0);

            T g_m_d = gamma + dhx*(cf0+cf1);

            T res = rhs(i,0,0) - gamma*phi(i,0,0)
                - dhx*(phi(i-1,0,0) + phi(i+1,0,0));

            phi(i,0,0) = phi(i,0,0) + res /g_m_d;
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_gsrb_m (int i, int, int, Array4<T> const& phi, Array4<T const> const& rhs,
                       T dhx,
                       Array4<T const> const& f0, Array4<int const> const& m0,
                       Array4<T const> const& f1, Array4<int const> const& m1,
                       Box const& vbox, int redblack, T dx, T probxlo) noexcept
{
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    if ((i+redblack)%2 == 0) {
        T cf0 = (i == vlo.x && m0(vlo.x-1,0,0) > 0)
            ? f0(vlo.x,0,0) : T(0.0);
        T cf1 = (i == vhi.x && m1(vhi.x+1,0,0) > 0)
            ? f1(vhi.x,0,0) : T(0.0);

        T rel = (probxlo + i   *dx) * (probxlo + i   *dx);
        T rer = (probxlo +(i+1)*dx) * (probxlo +(i+1)*dx);

        T gamma = -dhx*(rel+rer);

        T g_m_d = gamma + dhx*(rel*cf0+rer*cf1);

        T res = rhs(i,0,0) - gamma*phi(i,0,0)
            - dhx*(rel*phi(i-1,0,0) + rer*phi(i+1,0,0));

        phi(i,0,0) = phi(i,0,0) + res /g_m_d;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_normalize (int i, int, int, Array4<T> const& x,
                          T dhx, T dx, T probxlo) noexcept
{
    T rel = (probxlo + i   *dx) * (probxlo + i   *dx);
    T rer = (probxlo +(i+1)*dx) * (probxlo +(i+1)*dx);
    x(i,0,0) /= (-dhx*(rel+rer));
}

}

#endif
