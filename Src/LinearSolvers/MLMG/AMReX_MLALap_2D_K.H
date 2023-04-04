#ifndef AMREX_MLALAP_2D_K_H_
#define AMREX_MLALAP_2D_K_H_
#include <AMReX_Config.H>

namespace amrex::TwoD {

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_adotx (Box const& box, Array4<RT> const& y,
                   Array4<RT const> const& x,
                   Array4<RT const> const& a,
                   GpuArray<RT,2> const& dxinv,
                   RT alpha, RT beta, int ncomp) noexcept
{
    const RT dhx = beta*dxinv[0]*dxinv[0];
    const RT dhy = beta*dxinv[1]*dxinv[1];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                y(i,j,0,n) = alpha*a(i,j,0,n)*x(i,j,0,n)
                    - dhx * (x(i-1,j,0,n) - RT(2.)*x(i,j,0,n) + x(i+1,j,0,n))
                    - dhy * (x(i,j-1,0,n) - RT(2.)*x(i,j,0,n) + x(i,j+1,0,n));
            }
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_adotx_m (Box const& box, Array4<RT> const& y,
                     Array4<RT const> const& x,
                     Array4<RT const> const& a,
                     GpuArray<RT,2> const& dxinv,
                     RT alpha, RT beta, RT dx, RT probxlo, int ncomp) noexcept
{
    const RT dhx = beta*dxinv[0]*dxinv[0];
    const RT dhy = beta*dxinv[1]*dxinv[1];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                RT rel = probxlo + i*dx;
                RT rer = probxlo +(i+1)*dx;
                RT rc = probxlo + (i+RT(0.5))*dx;
                y(i,j,0,n) = alpha*a(i,j,0,n)*x(i,j,0,n)*rc
                    - dhx * (rer*(x(i+1,j,0,n) - x(i  ,j,0,n))
                           - rel*(x(i  ,j,0,n) - x(i-1,j,0,n)))
                    - dhy * rc * (x(i,j-1,0,n) - RT(2.)*x(i,j,0,n) + x(i,j+1,0,n));
            }
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_normalize (Box const& box, Array4<RT> const& x,
                       Array4<RT const> const& a,
                       GpuArray<RT,2> const& dxinv,
                       RT alpha, RT beta, int ncomp) noexcept
{
    const RT dhx = beta*dxinv[0]*dxinv[0];
    const RT dhy = beta*dxinv[1]*dxinv[1];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                x(i,j,0,n) /= alpha*a(i,j,0,n) + RT(2.0)*(dhx + dhy);
            }
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_normalize_m (Box const& box, Array4<RT> const& x,
                         Array4<RT const> const& a,
                         GpuArray<RT,2> const& dxinv,
                         RT alpha, RT beta,
                         RT dx, RT probxlo, int ncomp) noexcept
{
    const RT dhx = beta*dxinv[0]*dxinv[0];
    const RT dhy = beta*dxinv[1]*dxinv[1];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                RT rel = probxlo + i*dx;
                RT rer = probxlo +(i+1)*dx;
                RT rc = probxlo + (i+RT(0.5))*dx;
                x(i,j,0,n) /= alpha*a(i,j,0,n)*rc
                    + dhx*(rel+rer) + dhy*(rc*RT(2.0));
            }
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
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fx(i,j,0,n) = -fac*(sol(i,j,0,n)-sol(i-1,j,0,n));
            }
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
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                RT re = probxlo + i*dx;
                fx(i,j,0,n) = -fac*re*(sol(i,j,0,n)-sol(i-1,j,0,n));
            }
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_xface (Box const& box, Array4<RT> const& fx,
                        Array4<RT const> const& sol, RT fac, int xlen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            int i = lo.x;
            fx(i,j,0,n) = -fac*(sol(i,j,0,n)-sol(i-1,j,0,n));
            i += xlen;
            fx(i,j,0,n) = -fac*(sol(i,j,0,n)-sol(i-1,j,0,n));
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_xface_m (Box const& box, Array4<RT> const& fx,
                          Array4<RT const> const& sol, RT fac, int xlen,
                          RT dx, RT probxlo, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            int i = lo.x;
            RT re = probxlo + i*dx;
            fx(i,j,0,n) = -fac*re*(sol(i,j,0,n)-sol(i-1,j,0,n));
            i += xlen;
            re = probxlo + i*dx;
            fx(i,j,0,n) = -fac*re*(sol(i,j,0,n)-sol(i-1,j,0,n));
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_y (Box const& box, Array4<RT> const& fy,
                    Array4<RT const> const& sol, RT fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fy(i,j,0,n) = -fac*(sol(i,j,0,n)-sol(i,j-1,0,n));
            }
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_y_m (Box const& box, Array4<RT> const& fy,
                      Array4<RT const> const& sol, RT fac,
                      RT dx, RT probxlo, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                RT rc = probxlo + (i+RT(0.5))*dx;
                fy(i,j,0,n) = -fac*rc*(sol(i,j,0,n)-sol(i,j-1,0,n));
            }
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_yface (Box const& box, Array4<RT> const& fy,
                        Array4<RT const> const& sol, RT fac, int ylen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        int j = lo.y;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,0,n) = -fac*(sol(i,j,0,n)-sol(i,j-1,0,n));
        }
        j += ylen;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,0,n) = -fac*(sol(i,j,0,n)-sol(i,j-1,0,n));
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_yface_m (Box const& box, Array4<RT> const& fy,
                          Array4<RT const> const& sol, RT fac, int ylen,
                          RT dx, RT probxlo, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
        int j = lo.y;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            RT rc = probxlo + (i+RT(0.5))*dx;
            fy(i,j,0,n) = -fac*rc*(sol(i,j,0,n)-sol(i,j-1,0,n));
        }
        j += ylen;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            RT rc = probxlo + (i+RT(0.5))*dx;
            fy(i,j,0,n) = -fac*rc*(sol(i,j,0,n)-sol(i,j-1,0,n));
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_gsrb (Box const& box, Array4<RT> const& phi,
                  Array4<RT const> const& rhs, RT alpha,
                  RT dhx, RT dhy, Array4<RT const> const& a,
                  Array4<RT const> const& f0, Array4<int const> const& m0,
                  Array4<RT const> const& f1, Array4<int const> const& m1,
                  Array4<RT const> const& f2, Array4<int const> const& m2,
                  Array4<RT const> const& f3, Array4<int const> const& m3,
                  Box const& vbox, int redblack, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if ((i+j+redblack)%2 == 0) {
                    RT cf0 = (i == vlo.x && m0(vlo.x-1,j,0) > 0)
                        ? f0(vlo.x,j,0,n) : RT(0.0);
                    RT cf1 = (j == vlo.y && m1(i,vlo.y-1,0) > 0)
                        ? f1(i,vlo.y,0,n) : RT(0.0);
                    RT cf2 = (i == vhi.x && m2(vhi.x+1,j,0) > 0)
                        ? f2(vhi.x,j,0,n) : RT(0.0);
                    RT cf3 = (j == vhi.y && m3(i,vhi.y+1,0) > 0)
                        ? f3(i,vhi.y,0,n) : RT(0.0);

                    RT delta = dhx*(cf0 + cf2) + dhy*(cf1 + cf3);

                    RT gamma = alpha*a(i,j,0,n) + RT(2.0)*(dhx + dhy);

                    RT rho = dhx*(phi(i-1,j,0,n) + phi(i+1,j,0,n))
                        +    dhy*(phi(i,j-1,0,n) + phi(i,j+1,0,n));

                    phi(i,j,0,n) = (rhs(i,j,0,n) + rho - phi(i,j,0,n)*delta)
                        / (gamma - delta);
                }
            }
        }
    }
}

template <typename RT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_gsrb_m (Box const& box, Array4<RT> const& phi,
                    Array4<RT const> const& rhs, RT alpha,
                    RT dhx, RT dhy, Array4<RT const> const& a,
                    Array4<RT const> const& f0, Array4<int const> const& m0,
                    Array4<RT const> const& f1, Array4<int const> const& m1,
                    Array4<RT const> const& f2, Array4<int const> const& m2,
                    Array4<RT const> const& f3, Array4<int const> const& m3,
                    Box const& vbox, int redblack,
                    RT dx, RT probxlo, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if ((i+j+redblack)%2 == 0) {
                    RT cf0 = (i == vlo.x && m0(vlo.x-1,j,0) > 0)
                        ? f0(vlo.x,j,0,n) : RT(0.0);
                    RT cf1 = (j == vlo.y && m1(i,vlo.y-1,0) > 0)
                        ? f1(i,vlo.y,0,n) : RT(0.0);
                    RT cf2 = (i == vhi.x && m2(vhi.x+1,j,0) > 0)
                        ? f2(vhi.x,j,0,n) : RT(0.0);
                    RT cf3 = (j == vhi.y && m3(i,vhi.y+1,0) > 0)
                        ? f3(i,vhi.y,0,n) : RT(0.0);

                    RT rel = probxlo + i*dx;
                    RT rer = probxlo +(i+1)*dx;
                    RT rc = probxlo + (i+RT(0.5))*dx;

                    RT delta = dhx*(rel*cf0 + rer*cf2) + dhy*rc*(cf1 + cf3);

                    RT gamma = alpha*a(i,j,0,n)*rc + dhx*(rel+rer) + dhy*(RT(2.)*rc);

                    RT rho = dhx*(rel*phi(i-1,j,0,n) + rer*phi(i+1,j,0,n))
                      +      dhy* rc*(phi(i,j-1,0,n) +     phi(i,j+1,0,n));

                    phi(i,j,0,n) = (rhs(i,j,0,n) + rho - phi(i,j,0,n)*delta)
                        / (gamma - delta);
                }
            }
        }
    }
}

}

#if (AMREX_SPACEDIM == 2)
namespace amrex {
    using namespace TwoD;
}
#endif

#endif
