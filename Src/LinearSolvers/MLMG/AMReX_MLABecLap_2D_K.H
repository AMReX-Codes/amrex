#ifndef AMREX_MLABECLAP_2D_K_H_
#define AMREX_MLABECLAP_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_adotx (Box const& box, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Array4<Real const> const& a,
                      Array4<Real const> const& bX,
                      Array4<Real const> const& bY,
                      GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                      Real alpha, Real beta, int ncomp) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            y(i,j,0,n) = alpha*a(i,j,0)*x(i,j,0,n)
                - dhx * (bX(i+1,j,0)*(x(i+1,j,0,n) - x(i  ,j,0,n))
                       - bX(i  ,j,0)*(x(i  ,j,0,n) - x(i-1,j,0,n)))
                - dhy * (bY(i,j+1,0)*(x(i,j+1,0,n) - x(i,j  ,0,n))
                       - bY(i,j  ,0)*(x(i,j  ,0,n) - x(i,j-1,0,n)));
        }
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_normalize (Box const& box, Array4<Real> const& x,
                          Array4<Real const> const& a,
                          Array4<Real const> const& bX,
                          Array4<Real const> const& bY,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                          Real alpha, Real beta, int ncomp) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            x(i,j,0,n) /= alpha*a(i,j,0)
                + dhx*(bX(i,j,0)+bX(i+1,j,0))
                + dhy*(bY(i,j,0)+bY(i,j+1,0));
        }
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_flux_x (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                       Array4<Real const> const& bx, Real fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fx(i,j,0,n) = -fac*bx(i,j,0)*(sol(i,j,0,n)-sol(i-1,j,0,n));
        }
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_flux_xface (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                           Array4<Real const> const& bx, Real fac, int xlen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for     (int j = lo.y; j <= hi.y; ++j) {
        int i = lo.x;
        fx(i,j,0,n) = -fac*bx(i,j,0)*(sol(i,j,0,n)-sol(i-1,j,0,n));
        i += xlen;
        fx(i,j,0,n) = -fac*bx(i,j,0)*(sol(i,j,0,n)-sol(i-1,j,0,n));
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_flux_y (Box const& box, Array4<Real> const& fy, Array4<Real const> const& sol,
                       Array4<Real const> const& by, Real fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,0,n) = -fac*by(i,j,0)*(sol(i,j,0,n)-sol(i,j-1,0,n));
        }
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_flux_yface (Box const& box, Array4<Real> const& fy, Array4<Real const> const& sol,
                           Array4<Real const> const& by, Real fac, int ylen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    int j = lo.y;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fy(i,j,0,n) = -fac*by(i,j,0)*(sol(i,j,0,n)-sol(i,j-1,0,n));
    }
    j += ylen;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fy(i,j,0,n) = -fac*by(i,j,0)*(sol(i,j,0,n)-sol(i,j-1,0,n));
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void abec_gsrb (Box const& box, Array4<Real> const& phi,
                Array4<Real const> const& rhs, Real alpha,
                Real dhx, Real dhy, Array4<Real const> const& a,
                Array4<Real const> const& bX,
                Array4<Real const> const& bY,
                Array4<Real const> const& f0, Array4<int const> const& m0,
                Array4<Real const> const& f1, Array4<int const> const& m1,
                Array4<Real const> const& f2, Array4<int const> const& m2,
                Array4<Real const> const& f3, Array4<int const> const& m3,
                Box const& vbox, int nc, int redblack) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    for (int n = 0; n < nc; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if ((i+j+redblack)%2 == 0) {
                    Real cf0 = (i == vlo.x and m0(vlo.x-1,j,0) > 0)
                        ? f0(vlo.x,j,0,n) : 0.0;
                    Real cf1 = (j == vlo.y and m1(i,vlo.y-1,0) > 0)
                        ? f1(i,vlo.y,0,n) : 0.0;
                    Real cf2 = (i == vhi.x and m2(vhi.x+1,j,0) > 0)
                        ? f2(vhi.x,j,0,n) : 0.0;
                    Real cf3 = (j == vhi.y and m3(i,vhi.y+1,0) > 0)
                        ? f3(i,vhi.y,0,n) : 0.0;

                    Real delta = dhx*(bX(i,j,0)*cf0 + bX(i+1,j,0)*cf2)
                              +  dhy*(bY(i,j,0)*cf1 + bY(i,j+1,0)*cf3);

                    Real gamma = alpha*a(i,j,0)
                        +   dhx*( bX(i,j,0) + bX(i+1,j,0) )
                        +   dhy*( bY(i,j,0) + bY(i,j+1,0) );

                    Real rho = dhx*(bX(i  ,j  ,0)*phi(i-1,j  ,0,n)
                                  + bX(i+1,j  ,0)*phi(i+1,j  ,0,n))
                              +dhy*(bY(i  ,j  ,0)*phi(i  ,j-1,0,n)
                                  + bY(i  ,j+1,0)*phi(i  ,j+1,0,n));

                    phi(i,j,0,n) = (rhs(i,j,0,n) + rho - phi(i,j,0,n)*delta)
                        / (gamma - delta);
                }
            }
        }
    }
}

}
#endif
