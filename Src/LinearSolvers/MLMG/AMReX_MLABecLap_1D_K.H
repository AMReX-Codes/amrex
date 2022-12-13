#ifndef AMREX_MLABECLAP_1D_K_H_
#define AMREX_MLABECLAP_1D_K_H_
#include <AMReX_Config.H>

namespace amrex {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_adotx (int i, int, int, int n, Array4<T> const& y,
                      Array4<T const> const& x,
                      Array4<T const> const& a,
                      Array4<T const> const& bX,
                      GpuArray<T,AMREX_SPACEDIM> const& dxinv,
                      T alpha, T beta) noexcept
{
    const T dhx = beta*dxinv[0]*dxinv[0];
    y(i,0,0,n) = alpha*a(i,0,0)*x(i,0,0,n)
        - dhx * (bX(i+1,0,0,n)*(x(i+1,0,0,n) - x(i  ,0,0,n))
               - bX(i  ,0,0,n)*(x(i  ,0,0,n) - x(i-1,0,0,n)));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_adotx_os (int i, int, int, int n, Array4<T> const& y,
                         Array4<T const> const& x,
                         Array4<T const> const& a,
                         Array4<T const> const& bX,
                         Array4<int const> const& osm,
                         GpuArray<T,AMREX_SPACEDIM> const& dxinv,
                         T alpha, T beta) noexcept
{
    const T dhx = beta*dxinv[0]*dxinv[0];
    if (osm(i,0,0) == 0) {
        y(i,0,0,n) = T(0.0);
    } else {
        y(i,0,0,n) = alpha*a(i,0,0)*x(i,0,0,n)
            - dhx * (bX(i+1,0,0,n)*(x(i+1,0,0,n) - x(i  ,0,0,n))
                   - bX(i  ,0,0,n)*(x(i  ,0,0,n) - x(i-1,0,0,n)));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_normalize (int i, int, int, int n, Array4<T> const& x,
                          Array4<T const> const& a,
                          Array4<T const> const& bX,
                          GpuArray<T,AMREX_SPACEDIM> const& dxinv,
                          T alpha, T beta) noexcept
{
    const T dhx = beta*dxinv[0]*dxinv[0];
    x(i,0,0,n) /= alpha*a(i,0,0) + dhx*(bX(i,0,0,n)+bX(i+1,0,0,n));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_x (Box const& box, Array4<T> const& fx, Array4<T const> const& sol,
                       Array4<T const> const& bx, T fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fx(i,0,0,n) = -fac*bx(i,0,0,n)*(sol(i,0,0,n)-sol(i-1,0,0,n));
    }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_xface (Box const& box, Array4<T> const& fx, Array4<T const> const& sol,
                           Array4<T const> const& bx, T fac, int xlen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);

    for (int n = 0; n < ncomp; ++n) {
    int i = lo.x;
    fx(i,0,0,n) = -fac*bx(i,0,0,n)*(sol(i,0,0,n)-sol(i-1,0,0,n));
    i += xlen;
    fx(i,0,0,n) = -fac*bx(i,0,0,n)*(sol(i,0,0,n)-sol(i-1,0,0,n));
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void abec_gsrb (int i, int, int, int n, Array4<T> const& phi, Array4<T const> const& rhs,
                T alpha, Array4<T const> const& a,
                T dhx,
                Array4<T const> const& bX,
                Array4<int const> const& m0,
                Array4<int const> const& m1,
                Array4<T const> const& f0,
                Array4<T const> const& f1,
                Box const& vbox, int redblack) noexcept
{
    if ((i+redblack)%2 == 0) {
        const auto vlo = amrex::lbound(vbox);
        const auto vhi = amrex::ubound(vbox);

        T cf0 = (i == vlo.x && m0(vlo.x-1,0,0) > 0)
            ? f0(vlo.x,0,0,n) : T(0.0);
        T cf1 = (i == vhi.x && m1(vhi.x+1,0,0) > 0)
            ? f1(vhi.x,0,0,n) : T(0.0);

        T delta = dhx*(bX(i,0,0,n)*cf0 + bX(i+1,0,0,n)*cf1);

        T gamma = alpha*a(i,0,0)
            +   dhx*( bX(i,0,0,n) + bX(i+1,0,0,n) );

        T rho = dhx*(bX(i  ,0  ,0,n)*phi(i-1,0  ,0,n)
                        + bX(i+1,0  ,0,n)*phi(i+1,0  ,0,n));

        phi(i,0,0,n) = (rhs(i,0,0,n) + rho - phi(i,0,0,n)*delta)
            / (gamma - delta);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void abec_gsrb_os (int i, int, int, int n, Array4<T> const& phi, Array4<T const> const& rhs,
                   T alpha, Array4<T const> const& a,
                   T dhx,
                   Array4<T const> const& bX,
                   Array4<int const> const& m0,
                   Array4<int const> const& m1,
                   Array4<T const> const& f0,
                   Array4<T const> const& f1,
                   Array4<int const> const& osm,
                   Box const& vbox, int redblack) noexcept
{
    if ((i+redblack)%2 == 0) {
        if (osm(i,0,0) == 0) {
            phi(i,0,0) = T(0.0);
        } else {
            const auto vlo = amrex::lbound(vbox);
            const auto vhi = amrex::ubound(vbox);

            T cf0 = (i == vlo.x && m0(vlo.x-1,0,0) > 0)
                ? f0(vlo.x,0,0,n) : T(0.0);
            T cf1 = (i == vhi.x && m1(vhi.x+1,0,0) > 0)
                ? f1(vhi.x,0,0,n) : T(0.0);

            T delta = dhx*(bX(i,0,0,n)*cf0 + bX(i+1,0,0,n)*cf1);

            T gamma = alpha*a(i,0,0)
                +   dhx*( bX(i,0,0,n) + bX(i+1,0,0,n) );

            T rho = dhx*(bX(i  ,0  ,0,n)*phi(i-1,0  ,0,n)
                            + bX(i+1,0  ,0,n)*phi(i+1,0  ,0,n));

            phi(i,0,0,n) = (rhs(i,0,0,n) + rho - phi(i,0,0,n)*delta)
                / (gamma - delta);
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void abec_gsrb_with_line_solve (
                Box const& /*box*/, Array4<T> const& /*phi*/, Array4<T const> const& /*rhs*/,
                T /*alpha*/, Array4<T const> const& /*a*/,
                T /*dhx*/,
                Array4<T const> const& /*bX*/,
                Array4<int const> const& /*m0*/,
                Array4<int const> const& /*m1*/,
                Array4<T const> const& /*f0*/,
                Array4<T const> const& /*f1*/,
                Box const& /*vbox*/, int /*redblack*/, int /*nc*/) noexcept
{
    amrex::Abort("abec_gsrb_with_line_solve not implemented in 1D");
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void overset_rescale_bcoef_x (Box const& box, Array4<T> const& bX, Array4<int const> const& osm,
                              int ncomp, T osfac) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int n = 0; n < ncomp; ++n) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((osm(i-1,0,0)+osm(i,0,0)) == 1) {
                bX(i,0,0,n) *= osfac;
            }
        }
    }
}

}
#endif
