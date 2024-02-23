#ifndef AMREX_MLPOISSON_3D_K_H_
#define AMREX_MLPOISSON_3D_K_H_
#include <AMReX_Config.H>

namespace amrex {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx (int i, int j, int k, Array4<T> const& y,
                      Array4<T const> const& x,
                      T dhx, T dhy, T dhz) noexcept
{
    y(i,j,k) = dhx * (x(i-1,j,k) - T(2.0)*x(i,j,k) + x(i+1,j,k))
        +      dhy * (x(i,j-1,k) - T(2.0)*x(i,j,k) + x(i,j+1,k))
        +      dhz * (x(i,j,k-1) - T(2.0)*x(i,j,k) + x(i,j,k+1));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx_os (int i, int j, int k, Array4<T> const& y,
                         Array4<T const> const& x,
                         Array4<int const> const& osm,
                         T dhx, T dhy, T dhz) noexcept
{
    if (osm(i,j,k) == 0) {
        y(i,j,k) = T(0.0);
    } else {
        y(i,j,k) = dhx * (x(i-1,j,k) - T(2.0)*x(i,j,k) + x(i+1,j,k))
            +      dhy * (x(i,j-1,k) - T(2.0)*x(i,j,k) + x(i,j+1,k))
            +      dhz * (x(i,j,k-1) - T(2.0)*x(i,j,k) + x(i,j,k+1));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_x (Box const& box, Array4<T> const& fx,
                       Array4<T const> const& sol, T dxinv) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fx(i,j,k) = dxinv*(sol(i,j,k)-sol(i-1,j,k));
            }
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

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            int i = lo.x;
            fx(i,j,k) = dxinv*(sol(i,j,k)-sol(i-1,j,k));
            i += xlen;
            fx(i,j,k) = dxinv*(sol(i,j,k)-sol(i-1,j,k));
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_y (Box const& box, Array4<T> const& fy,
                       Array4<T const> const& sol, T dyinv) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fy(i,j,k) = dyinv*(sol(i,j,k)-sol(i,j-1,k));
            }
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

    for     (int k = lo.z; k <= hi.z; ++k) {
        int j = lo.y;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,k) = dyinv*(sol(i,j,k)-sol(i,j-1,k));
        }
        j += ylen;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,k) = dyinv*(sol(i,j,k)-sol(i,j-1,k));
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_z (Box const& box, Array4<T> const& fz,
                       Array4<T const> const& sol, T dzinv) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fz(i,j,k) = dzinv*(sol(i,j,k)-sol(i,j,k-1));
            }
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_zface (Box const& box, Array4<T> const& fz,
                           Array4<T const> const& sol, T dzinv, int zlen) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = lo.z;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fz(i,j,k) = dzinv*(sol(i,j,k)-sol(i,j,k-1));
        }
    }

    k += zlen;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fz(i,j,k) = dzinv*(sol(i,j,k)-sol(i,j,k-1));
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_gsrb (int i, int j, int k, Array4<T> const& phi,
                     Array4<T const> const& rhs,
                     T dhx, T dhy, T dhz,
                     Array4<T const> const& f0, Array4<int const> const& m0,
                     Array4<T const> const& f1, Array4<int const> const& m1,
                     Array4<T const> const& f2, Array4<int const> const& m2,
                     Array4<T const> const& f3, Array4<int const> const& m3,
                     Array4<T const> const& f4, Array4<int const> const& m4,
                     Array4<T const> const& f5, Array4<int const> const& m5,
                     Box const& vbox, int redblack) noexcept
{
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    constexpr T omega = T(1.15);

    const T gamma = T(-2.)*(dhx+dhy+dhz);

    if ((i+j+k+redblack)%2 == 0) {
        T cf0 = (i == vlo.x && m0(vlo.x-1,j,k) > 0)
            ? f0(vlo.x,j,k) : T(0.0);
        T cf1 = (j == vlo.y && m1(i,vlo.y-1,k) > 0)
            ? f1(i,vlo.y,k) : T(0.0);
        T cf2 = (k == vlo.z && m2(i,j,vlo.z-1) > 0)
            ? f2(i,j,vlo.z) : T(0.0);
        T cf3 = (i == vhi.x && m3(vhi.x+1,j,k) > 0)
            ? f3(vhi.x,j,k) : T(0.0);
        T cf4 = (j == vhi.y && m4(i,vhi.y+1,k) > 0)
            ? f4(i,vhi.y,k) : T(0.0);
        T cf5 = (k == vhi.z && m5(i,j,vhi.z+1) > 0)
            ? f5(i,j,vhi.z) : T(0.0);

        T g_m_d = gamma + dhx*(cf0+cf3) + dhy*(cf1+cf4) + dhz*(cf2+cf5);

        T res = rhs(i,j,k) - gamma*phi(i,j,k)
            - dhx*(phi(i-1,j,k) + phi(i+1,j,k))
            - dhy*(phi(i,j-1,k) + phi(i,j+1,k))
            - dhz*(phi(i,j,k-1) + phi(i,j,k+1));

        phi(i,j,k) = phi(i,j,k) + omega/g_m_d * res;
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_gsrb_os (int i, int j, int k, Array4<T> const& phi,
                        Array4<T const> const& rhs,
                        Array4<int const> const& osm,
                        T dhx, T dhy, T dhz,
                        Array4<T const> const& f0, Array4<int const> const& m0,
                        Array4<T const> const& f1, Array4<int const> const& m1,
                        Array4<T const> const& f2, Array4<int const> const& m2,
                        Array4<T const> const& f3, Array4<int const> const& m3,
                        Array4<T const> const& f4, Array4<int const> const& m4,
                        Array4<T const> const& f5, Array4<int const> const& m5,
                        Box const& vbox, int redblack) noexcept
{
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    constexpr T omega = T(1.15);

    const T gamma = T(-2.)*(dhx+dhy+dhz);

    if ((i+j+k+redblack)%2 == 0) {
        if (osm(i,j,k) == 0) {
            phi(i,j,k) = T(0.0);
        } else {
            T cf0 = (i == vlo.x && m0(vlo.x-1,j,k) > 0)
                ? f0(vlo.x,j,k) : T(0.0);
            T cf1 = (j == vlo.y && m1(i,vlo.y-1,k) > 0)
                ? f1(i,vlo.y,k) : T(0.0);
            T cf2 = (k == vlo.z && m2(i,j,vlo.z-1) > 0)
                ? f2(i,j,vlo.z) : T(0.0);
            T cf3 = (i == vhi.x && m3(vhi.x+1,j,k) > 0)
                ? f3(vhi.x,j,k) : T(0.0);
            T cf4 = (j == vhi.y && m4(i,vhi.y+1,k) > 0)
                ? f4(i,vhi.y,k) : T(0.0);
            T cf5 = (k == vhi.z && m5(i,j,vhi.z+1) > 0)
                ? f5(i,j,vhi.z) : T(0.0);

            T g_m_d = gamma + dhx*(cf0+cf3) + dhy*(cf1+cf4) + dhz*(cf2+cf5);

            T res = rhs(i,j,k) - gamma*phi(i,j,k)
                - dhx*(phi(i-1,j,k) + phi(i+1,j,k))
                - dhy*(phi(i,j-1,k) + phi(i,j+1,k))
                - dhz*(phi(i,j,k-1) + phi(i,j,k+1));

            phi(i,j,k) = phi(i,j,k) + omega/g_m_d * res;
        }
    }
}

}

#endif
