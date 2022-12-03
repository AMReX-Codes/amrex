#ifndef AMREX_MLALAP_1D_K_H_
#define AMREX_MLALAP_1D_K_H_
#include <AMReX_Config.H>

namespace amrex {

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_adotx (Box const& box, Array4<RT> const& y,
                   Array4<RT const> const& x,
                   Array4<RT const> const& a,
                   GpuArray<RT,AMREX_SPACEDIM> const& dxinv,
                   RT alpha, RT beta, int ncomp) noexcept
{
    const RT dhx = beta*dxinv[0]*dxinv[0];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            y(i,0,0,n) = alpha*a(i,0,0,n)*x(i,0,0,n)
                - dhx * (x(i+1,0,0,n) - RT(2.)*x(i,0,0,n) + x(i-1,0,0,n));
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_adotx_m (Box const& box, Array4<RT> const& y,
                     Array4<RT const> const& x,
                     Array4<RT const> const& a,
                     GpuArray<RT,AMREX_SPACEDIM> const& dxinv,
                     RT alpha, RT beta,
                     RT dx, RT probxlo, int ncomp) noexcept
{
    const RT dhx = beta*dxinv[0]*dxinv[0];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            RT rel = (probxlo + i   *dx) * (probxlo + i   *dx);
            RT rer = (probxlo +(i+1)*dx) * (probxlo +(i+1)*dx);
            RT rc  = (probxlo +(i+RT(.5))*dx) * (probxlo +(i+RT(.5))*dx);
            y(i,0,0,n) = alpha*a(i,0,0,n)*x(i,0,0,n)*rc
                - dhx * (rer*(x(i+1,0,0,n) - x(i  ,0,0,n))
                       - rel*(x(i  ,0,0,n) - x(i-1,0,0,n)));
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_normalize (Box const& box, Array4<RT> const& x,
                       Array4<RT const> const& a,
                       GpuArray<RT,AMREX_SPACEDIM> const& dxinv,
                       RT alpha, RT beta, int ncomp) noexcept
{
    const RT dhx = beta*dxinv[0]*dxinv[0];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            x(i,0,0,n) /= alpha*a(i,0,0,n) + dhx*RT(2.0);
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_normalize_m (Box const& box, Array4<RT> const& x,
                         Array4<RT const> const& a,
                         GpuArray<RT,AMREX_SPACEDIM> const& dxinv,
                         RT alpha, RT beta, RT dx, RT probxlo, int ncomp) noexcept
{
    const RT dhx = beta*dxinv[0]*dxinv[0];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            RT rel = (probxlo + i   *dx) * (probxlo + i   *dx);
            RT rer = (probxlo +(i+1)*dx) * (probxlo +(i+1)*dx);
            RT rc  = (probxlo +(i+RT(.5))*dx) * (probxlo +(i+RT(.5))*dx);
            x(i,0,0,n) /= alpha*a(i,0,0,n)*rc + dhx*(rel+rer);
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_x (Box const& box, Array4<RT> const& fx,
                    Array4<RT const> const& sol, RT fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fx(i,0,0,n) = -fac*(sol(i,0,0,n)-sol(i-1,0,0,n));
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_x_m (Box const& box, Array4<RT> const& fx,
                      Array4<RT const> const& sol, RT fac,
                      RT dx, RT probxlo, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            RT re = (probxlo + i*dx) * (probxlo + i*dx);
            fx(i,0,0,n) = -fac*re*(sol(i,0,0,n)-sol(i-1,0,0,n));
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_xface (Box const& box, Array4<RT> const& fx,
                        Array4<RT const> const& sol, RT fac, int xlen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);

    for (int n = 0; n < ncomp; ++n) {
        int i = lo.x;
        fx(i,0,0,n) = -fac*(sol(i,0,0,n)-sol(i-1,0,0,n));
        i += xlen;
        fx(i,0,0,n) = -fac*(sol(i,0,0,n)-sol(i-1,0,0,n));
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_xface_m (Box const& box, Array4<RT> const& fx,
                          Array4<RT const> const& sol, RT fac, int xlen,
                          RT dx, RT probxlo, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);

    for (int n = 0; n < ncomp; ++n) {
        int i = lo.x;
        RT re = (probxlo + i*dx) * (probxlo + i*dx);
        fx(i,0,0,n) = -fac*re*(sol(i,0,0,n)-sol(i-1,0,0,n));
        i += xlen;
        re = (probxlo + i*dx) * (probxlo + i*dx);
        fx(i,0,0,n) = -fac*re*(sol(i,0,0,n)-sol(i-1,0,0,n));
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_gsrb (Box const& box, Array4<RT> const& phi,
                  Array4<RT const> const& rhs, RT alpha,
                  RT dhx, Array4<RT const> const& a,
                  Array4<RT const> const& f0, Array4<int const> const& m0,
                  Array4<RT const> const& f1, Array4<int const> const& m1,
                  Box const& vbox, int redblack, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((i+redblack)%2 == 0) {
                RT cf0 = (i == vlo.x && m0(vlo.x-1,0,0) > 0)
                    ? f0(vlo.x,0,0,n) : RT(0.0);
                RT cf1 = (i == vhi.x && m1(vhi.x+1,0,0) > 0)
                    ? f1(vhi.x,0,0,n) : RT(0.0);

                RT delta = dhx*(cf0+cf1);

                RT gamma = alpha*a(i,0,0,n) + dhx*RT(2.0);

                RT rho = dhx*(phi(i-1,0,0,n) + phi(i+1,0,0,n));

                phi(i,0,0,n) = (rhs(i,0,0,n) + rho - phi(i,0,0,n)*delta)
                    / (gamma - delta);
            }
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_gsrb_m (Box const& box, Array4<RT> const& phi,
                    Array4<RT const> const& rhs, RT alpha,
                    RT dhx, Array4<RT const> const& a,
                    Array4<RT const> const& f0, Array4<int const> const& m0,
                    Array4<RT const> const& f1, Array4<int const> const& m1,
                    Box const& vbox, int redblack,
                    RT dx, RT probxlo, int ncomp)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((i+redblack)%2 == 0) {
                RT cf0 = (i == vlo.x && m0(vlo.x-1,0,0) > 0)
                    ? f0(vlo.x,0,0,n) : RT(0.0);
                RT cf1 = (i == vhi.x && m1(vhi.x+1,0,0) > 0)
                    ? f1(vhi.x,0,0,n) : RT(0.0);

                RT rel = (probxlo + i   *dx) * (probxlo + i   *dx);
                RT rer = (probxlo +(i+1)*dx) * (probxlo +(i+1)*dx);
                RT rc  = (probxlo +(i+RT(.5))*dx) * (probxlo +(i+RT(.5))*dx);

                RT delta = dhx*(rel*cf0 + rer*cf1);

                RT gamma = alpha*a(i,0,0,n)*rc + dhx*(rel+rer);

                RT rho = dhx*(rel*phi(i-1,0,0,n) + rer*phi(i+1,0,0,n));

                phi(i,0,0,n) = (rhs(i,0,0,n) + rho - phi(i,0,0,n)*delta)
                    / (gamma - delta);
            }
        }
    }
}

}
#endif
