#ifndef AMREX_MLPOISSON_2D_K_H_
#define AMREX_MLPOISSON_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx (int i, int j, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Real dhx, Real dhy) noexcept
{
    y(i,j,0) = dhx * (x(i-1,j,0) - 2.*x(i,j,0) + x(i+1,j,0))
        +      dhy * (x(i,j-1,0) - 2.*x(i,j,0) + x(i,j+1,0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx_m (int i, int j, Array4<Real> const& y,
                        Array4<Real const> const& x,
                        Real dhx, Real dhy, Real dx, Real probxlo) noexcept
{
    Real rel = probxlo + i*dx;
    Real rer = probxlo +(i+1)*dx;
    Real rc = probxlo + (i+0.5)*dx;
    y(i,j,0) = dhx * (rel*x(i-1,j,0) - (rel+rer)*x(i,j,0) + rer*x(i+1,j,0))
        +      dhy * rc *(x(i,j-1,0) -        2.*x(i,j,0) +     x(i,j+1,0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_x (Box const& box, Array4<Real> const& fx,
                       Array4<Real const> const& sol, Real dxinv) noexcept
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_x_m (Box const& box, Array4<Real> const& fx,
                         Array4<Real const> const& sol, Real dxinv,
                         Real dx, Real probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real re = probxlo + i*dx;
            fx(i,j,0) = dxinv*re*(sol(i,j,0)-sol(i-1,j,0));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_xface (Box const& box, Array4<Real> const& fx,
                           Array4<Real const> const& sol, Real dxinv, int xlen) noexcept
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_xface_m (Box const& box, Array4<Real> const& fx,
                             Array4<Real const> const& sol, Real dxinv, int xlen,
                             Real dx, Real probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        int i = lo.x;
        Real re = probxlo + i*dx;
        fx(i,j,0) = dxinv*re*(sol(i,j,0)-sol(i-1,j,0));
        i += xlen;
        re = probxlo + i*dx;
        fx(i,j,0) = dxinv*re*(sol(i,j,0)-sol(i-1,j,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_y (Box const& box, Array4<Real> const& fy,
                       Array4<Real const> const& sol, Real dyinv) noexcept
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_y_m (Box const& box, Array4<Real> const& fy,
                         Array4<Real const> const& sol, Real dyinv,
                         Real dx, Real probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real rc = probxlo + (i+0.5)*dx;
            fy(i,j,0) = dyinv*rc*(sol(i,j,0)-sol(i,j-1,0));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_yface (Box const& box, Array4<Real> const& fy,
                           Array4<Real const> const& sol, Real dyinv, int ylen) noexcept
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_flux_yface_m (Box const& box, Array4<Real> const& fy,
                             Array4<Real const> const& sol, Real dyinv, int ylen,
                             Real dx, Real probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int j = lo.y;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        Real rc = probxlo + (i+0.5)*dx;
        fy(i,j,0) = dyinv*rc*(sol(i,j,0)-sol(i,j-1,0));
    }
    j += ylen;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        Real rc = probxlo + (i+0.5)*dx;
        fy(i,j,0) = dyinv*rc*(sol(i,j,0)-sol(i,j-1,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_gsrb (Box const& box, Array4<Real> const& phi, Array4<Real const> const& rhs,
                     Real dhx, Real dhy,
                     Array4<Real const> const& f0, Array4<int const> const& m0,
                     Array4<Real const> const& f1, Array4<int const> const& m1,
                     Array4<Real const> const& f2, Array4<int const> const& m2,
                     Array4<Real const> const& f3, Array4<int const> const& m3,
                     Box const& vbox, int redblack) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    Real gamma = -2.0*dhx - 2.0*dhy;

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((i+j+redblack)%2 == 0) {
                Real cf0 = (i == vlo.x and m0(vlo.x-1,j,0) > 0)
                    ? f0(vlo.x,j,0) : 0.0;
                Real cf1 = (j == vlo.y and m1(i,vlo.y-1,0) > 0)
                    ? f1(i,vlo.y,0) : 0.0;
                Real cf2 = (i == vhi.x and m2(vhi.x+1,j,0) > 0)
                    ? f2(vhi.x,j,0) : 0.0;
                Real cf3 = (j == vhi.y and m3(i,vhi.y+1,0) > 0)
                    ? f3(i,vhi.y,0) : 0.0;

                Real g_m_d = gamma + dhx*(cf0+cf2) + dhy*(cf1+cf3);

                Real res = rhs(i,j,0) - gamma*phi(i,j,0)
                    - dhx*(phi(i-1,j,0) + phi(i+1,j,0))
                    - dhy*(phi(i,j-1,0) + phi(i,j+1,0));

                phi(i,j,0) = phi(i,j,0) + res /g_m_d;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_gsrb_m (Box const& box, Array4<Real> const& phi, Array4<Real const> const& rhs,
                       Real dhx, Real dhy,
                       Array4<Real const> const& f0, Array4<int const> const& m0,
                       Array4<Real const> const& f1, Array4<int const> const& m1,
                       Array4<Real const> const& f2, Array4<int const> const& m2,
                       Array4<Real const> const& f3, Array4<int const> const& m3,
                       Box const& vbox, int redblack, Real dx, Real probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((i+j+redblack)%2 == 0) {
                Real cf0 = (i == vlo.x and m0(vlo.x-1,j,0) > 0)
                    ? f0(vlo.x,j,0) : 0.0;
                Real cf1 = (j == vlo.y and m1(i,vlo.y-1,0) > 0)
                    ? f1(i,vlo.y,0) : 0.0;
                Real cf2 = (i == vhi.x and m2(vhi.x+1,j,0) > 0)
                    ? f2(vhi.x,j,0) : 0.0;
                Real cf3 = (j == vhi.y and m3(i,vhi.y+1,0) > 0)
                    ? f3(i,vhi.y,0) : 0.0;

                Real rel = probxlo + i*dx;
                Real rer = probxlo +(i+1)*dx;
                Real rc = probxlo + (i+0.5)*dx;

                Real gamma = -dhx*(rel+rer) - 2.0*dhy*rc;

                Real g_m_d = gamma + dhx*(rel*cf0+rer*cf2) + dhy*rc*(cf1+cf3);

                Real res = rhs(i,j,0) - gamma*phi(i,j,0)
                    - dhx*(rel*phi(i-1,j,0) + rer*phi(i+1,j,0))
                    - dhy*rc *(phi(i,j-1,0) +     phi(i,j+1,0));

                phi(i,j,0) = phi(i,j,0) + res /g_m_d;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_normalize (Box const& box, Array4<Real> const& x,
                          Real dhx, Real dhy, Real dx, Real probxlo) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real rel = probxlo + i*dx;
            Real rer = probxlo +(i+1)*dx;
            Real rc = probxlo + (i+0.5)*dx;
            x(i,j,0) /= (-dhx*(rel+rer) - dhy*rc*2.0);
        }
    }
}

}

#endif
