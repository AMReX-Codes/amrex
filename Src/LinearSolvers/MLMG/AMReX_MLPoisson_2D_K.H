#ifndef AMREX_MLPOISSON_2D_K_H_
#define AMREX_MLPOISSON_2D_K_H_
#include <AMReX_Config.H>

namespace amrex {

#if (AMREX_SPACEDIM != 2)
namespace TwoD {
#endif

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx (int i, int j, Array4<T> const& y,
                      Array4<T const> const& x,
                      T dhx, T dhy) noexcept
{
    y(i,j,0) = dhx * (x(i-1,j,0) - T(2.)*x(i,j,0) + x(i+1,j,0))
        +      dhy * (x(i,j-1,0) - T(2.)*x(i,j,0) + x(i,j+1,0));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx_os (int i, int j, Array4<T> const& y,
                         Array4<T const> const& x,
                         Array4<int const> const& osm,
                         T dhx, T dhy) noexcept
{
    if (osm(i,j,0) == 0) {
        y(i,j,0) = T(0.0);
    } else {
        y(i,j,0) = dhx * (x(i-1,j,0) - T(2.)*x(i,j,0) + x(i+1,j,0))
            +      dhy * (x(i,j-1,0) - T(2.)*x(i,j,0) + x(i,j+1,0));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx_m (int i, int j, Array4<T> const& y,
                        Array4<T const> const& x,
                        T dhx, T dhy, T dx, T probxlo) noexcept
{
    T rel = probxlo + i*dx;
    T rer = probxlo +(i+1)*dx;
    T rc = probxlo + (i+T(0.5))*dx;
    y(i,j,0) = dhx * (rel*x(i-1,j,0) - (rel+rer)*x(i,j,0) + rer*x(i+1,j,0))
        +      dhy * rc *(x(i,j-1,0) -  T(2.)*x(i,j,0) +     x(i,j+1,0));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_x (Box const& box, Array4<T> const& fx,
                       Array4<T const> const& sol, T dxinv) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fx(i,j,0) = dxinv*(sol(i,j,0)-sol(i-1,j,0));
        }
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

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            T re = probxlo + i*dx;
            fx(i,j,0) = dxinv*re*(sol(i,j,0)-sol(i-1,j,0));
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_xface (Box const& box, Array4<T> const& fx,
                           Array4<T const> const& sol, T dxinv, int xlen) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        int i = lo.x;
        fx(i,j,0) = dxinv*(sol(i,j,0)-sol(i-1,j,0));
        i += xlen;
        fx(i,j,0) = dxinv*(sol(i,j,0)-sol(i-1,j,0));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_xface_m (Box const& box, Array4<T> const& fx,
                             Array4<T const> const& sol, T dxinv, int xlen,
                             T dx, T probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        int i = lo.x;
        T re = probxlo + i*dx;
        fx(i,j,0) = dxinv*re*(sol(i,j,0)-sol(i-1,j,0));
        i += xlen;
        re = probxlo + i*dx;
        fx(i,j,0) = dxinv*re*(sol(i,j,0)-sol(i-1,j,0));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_y (Box const& box, Array4<T> const& fy,
                       Array4<T const> const& sol, T dyinv) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,0) = dyinv*(sol(i,j,0)-sol(i,j-1,0));
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_y_m (Box const& box, Array4<T> const& fy,
                         Array4<T const> const& sol, T dyinv,
                         T dx, T probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            T rc = probxlo + (i+T(0.5))*dx;
            fy(i,j,0) = dyinv*rc*(sol(i,j,0)-sol(i,j-1,0));
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_yface (Box const& box, Array4<T> const& fy,
                           Array4<T const> const& sol, T dyinv, int ylen) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int j = lo.y;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fy(i,j,0) = dyinv*(sol(i,j,0)-sol(i,j-1,0));
    }
    j += ylen;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fy(i,j,0) = dyinv*(sol(i,j,0)-sol(i,j-1,0));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_yface_m (Box const& box, Array4<T> const& fy,
                             Array4<T const> const& sol, T dyinv, int ylen,
                             T dx, T probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int j = lo.y;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        T rc = probxlo + (i+T(0.5))*dx;
        fy(i,j,0) = dyinv*rc*(sol(i,j,0)-sol(i,j-1,0));
    }
    j += ylen;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        T rc = probxlo + (i+T(0.5))*dx;
        fy(i,j,0) = dyinv*rc*(sol(i,j,0)-sol(i,j-1,0));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_gsrb (int i, int j, int, Array4<T> const& phi, Array4<T const> const& rhs,
                     T dhx, T dhy,
                     Array4<T const> const& f0, Array4<int const> const& m0,
                     Array4<T const> const& f1, Array4<int const> const& m1,
                     Array4<T const> const& f2, Array4<int const> const& m2,
                     Array4<T const> const& f3, Array4<int const> const& m3,
                     Box const& vbox, int redblack) noexcept
{
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    T gamma = T(-2.0)*(dhx+dhy);

    if ((i+j+redblack)%2 == 0) {
        T cf0 = (i == vlo.x && m0(vlo.x-1,j,0) > 0)
            ? f0(vlo.x,j,0) : T(0.0);
        T cf1 = (j == vlo.y && m1(i,vlo.y-1,0) > 0)
            ? f1(i,vlo.y,0) : T(0.0);
        T cf2 = (i == vhi.x && m2(vhi.x+1,j,0) > 0)
            ? f2(vhi.x,j,0) : T(0.0);
        T cf3 = (j == vhi.y && m3(i,vhi.y+1,0) > 0)
            ? f3(i,vhi.y,0) : T(0.0);

        T g_m_d = gamma + dhx*(cf0+cf2) + dhy*(cf1+cf3);

        T res = rhs(i,j,0) - gamma*phi(i,j,0)
            - dhx*(phi(i-1,j,0) + phi(i+1,j,0))
            - dhy*(phi(i,j-1,0) + phi(i,j+1,0));

        phi(i,j,0) = phi(i,j,0) + res /g_m_d;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_gsrb_os (int i, int j, int, Array4<T> const& phi, Array4<T const> const& rhs,
                        Array4<int const> const& osm, T dhx, T dhy,
                        Array4<T const> const& f0, Array4<int const> const& m0,
                        Array4<T const> const& f1, Array4<int const> const& m1,
                        Array4<T const> const& f2, Array4<int const> const& m2,
                        Array4<T const> const& f3, Array4<int const> const& m3,
                        Box const& vbox, int redblack) noexcept
{
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    T gamma = T(-2.0)*(dhx+dhy);

    if ((i+j+redblack)%2 == 0) {
        if (osm(i,j,0) == 0) {
            phi(i,j,0) = T(0.0);
        } else {
            T cf0 = (i == vlo.x && m0(vlo.x-1,j,0) > 0)
                ? f0(vlo.x,j,0) : T(0.0);
            T cf1 = (j == vlo.y && m1(i,vlo.y-1,0) > 0)
                ? f1(i,vlo.y,0) : T(0.0);
            T cf2 = (i == vhi.x && m2(vhi.x+1,j,0) > 0)
                ? f2(vhi.x,j,0) : T(0.0);
            T cf3 = (j == vhi.y && m3(i,vhi.y+1,0) > 0)
                ? f3(i,vhi.y,0) : T(0.0);

            T g_m_d = gamma + dhx*(cf0+cf2) + dhy*(cf1+cf3);

            T res = rhs(i,j,0) - gamma*phi(i,j,0)
                - dhx*(phi(i-1,j,0) + phi(i+1,j,0))
                - dhy*(phi(i,j-1,0) + phi(i,j+1,0));

            phi(i,j,0) = phi(i,j,0) + res /g_m_d;
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_gsrb_m (int i, int j, int, Array4<T> const& phi, Array4<T const> const& rhs,
                       T dhx, T dhy,
                       Array4<T const> const& f0, Array4<int const> const& m0,
                       Array4<T const> const& f1, Array4<int const> const& m1,
                       Array4<T const> const& f2, Array4<int const> const& m2,
                       Array4<T const> const& f3, Array4<int const> const& m3,
                       Box const& vbox, int redblack, T dx, T probxlo) noexcept
{
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    if ((i+j+redblack)%2 == 0) {
        T cf0 = (i == vlo.x && m0(vlo.x-1,j,0) > 0)
            ? f0(vlo.x,j,0) : T(0.0);
        T cf1 = (j == vlo.y && m1(i,vlo.y-1,0) > 0)
            ? f1(i,vlo.y,0) : T(0.0);
        T cf2 = (i == vhi.x && m2(vhi.x+1,j,0) > 0)
            ? f2(vhi.x,j,0) : T(0.0);
        T cf3 = (j == vhi.y && m3(i,vhi.y+1,0) > 0)
            ? f3(i,vhi.y,0) : T(0.0);

        T rel = probxlo + i*dx;
        T rer = probxlo +(i+1)*dx;
        T rc = probxlo + (i+T(0.5))*dx;

        T gamma = -dhx*(rel+rer) - T(2.0)*dhy*rc;

        T g_m_d = gamma + dhx*(rel*cf0+rer*cf2) + dhy*rc*(cf1+cf3);

        T res = rhs(i,j,0) - gamma*phi(i,j,0)
            - dhx*(rel*phi(i-1,j,0) + rer*phi(i+1,j,0))
            - dhy*rc *(phi(i,j-1,0) +     phi(i,j+1,0));

        phi(i,j,0) = phi(i,j,0) + res /g_m_d;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_normalize (int i, int j, int, Array4<T> const& x,
                          T dhx, T dhy, T dx, T probxlo) noexcept
{
    T rel = probxlo + i*dx;
    T rer = probxlo +(i+1)*dx;
    T rc = probxlo + (i+T(0.5))*dx;
    x(i,j,0) /= (-dhx*(rel+rer) - dhy*rc*T(2.0));
}

#if (AMREX_SPACEDIM != 2)
}
#endif

}

#endif
