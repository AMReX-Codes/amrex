#ifndef AMREX_MLABECLAP_3D_K_H_
#define AMREX_MLABECLAP_3D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_adotx (Box const& box, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Array4<Real const> const& a,
                      Array4<Real const> const& bX,
                      Array4<Real const> const& bY,
                      Array4<Real const> const& bZ,
                      GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                      Real alpha, Real beta)
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];
    const Real dhz = beta*dxinv[2]*dxinv[2];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                y(i,j,k) = alpha*a(i,j,k)*x(i,j,k)
                    - dhx * (bX(i+1,j,k)*(x(i+1,j,k) - x(i  ,j,k))
                           - bX(i  ,j,k)*(x(i  ,j,k) - x(i-1,j,k)))
                    - dhy * (bY(i,j+1,k)*(x(i,j+1,k) - x(i,j  ,k))
                           - bY(i,j  ,k)*(x(i,j  ,k) - x(i,j-1,k)))
                    - dhz * (bZ(i,j,k+1)*(x(i,j,k+1) - x(i,j,k  ))
                           - bZ(i,j,k  )*(x(i,j,k  ) - x(i,j,k-1)));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_normalize (Box const& box, Array4<Real> const& x,
                          Array4<Real const> const& a,
                          Array4<Real const> const& bX,
                          Array4<Real const> const& bY,
                          Array4<Real const> const& bZ,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                          Real alpha, Real beta)
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];
    const Real dhz = beta*dxinv[2]*dxinv[2];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                x(i,j,k) /= alpha*a(i,j,k)
                    + dhx*(bX(i,j,k)+bX(i+1,j,k))
                    + dhy*(bY(i,j,k)+bY(i,j+1,k))
                    + dhz*(bZ(i,j,k)+bZ(i,j,k+1));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_flux_x (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                       Array4<Real const> const& bx, Real fac)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fx(i,j,k) = -fac*bx(i,j,k)*(sol(i,j,k)-sol(i-1,j,k));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_flux_xface (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                           Array4<Real const> const& bx, Real fac, int xlen)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            int i = lo.x;
            fx(i,j,k) = -fac*bx(i,j,k)*(sol(i,j,k)-sol(i-1,j,k));
            i += xlen;
            fx(i,j,k) = -fac*bx(i,j,k)*(sol(i,j,k)-sol(i-1,j,k));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_flux_y (Box const& box, Array4<Real> const& fy, Array4<Real const> const& sol,
                       Array4<Real const> const& by, Real fac)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fy(i,j,k) = -fac*by(i,j,k)*(sol(i,j,k)-sol(i,j-1,k));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_flux_yface (Box const& box, Array4<Real> const& fy, Array4<Real const> const& sol,
                           Array4<Real const> const& by, Real fac, int ylen)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int k = lo.z; k <= hi.z; ++k) {
        int j = lo.y;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,k) = -fac*by(i,j,k)*(sol(i,j,k)-sol(i,j-1,k));
        }
        j += ylen;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,k) = -fac*by(i,j,k)*(sol(i,j,k)-sol(i,j-1,k));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_flux_z (Box const& box, Array4<Real> const& fz, Array4<Real const> const& sol,
                       Array4<Real const> const& bz, Real fac)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fz(i,j,k) = -fac*bz(i,j,k)*(sol(i,j,k)-sol(i,j,k-1));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_flux_zface (Box const& box, Array4<Real> const& fz, Array4<Real const> const& sol,
                           Array4<Real const> const& bz, Real fac, Real zlen)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = lo.z;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fz(i,j,k) = -fac*bz(i,j,k)*(sol(i,j,k)-sol(i,j,k-1));
        }
    }

    k += zlen;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fz(i,j,k) = -fac*bz(i,j,k)*(sol(i,j,k)-sol(i,j,k-1));
        }
    }
}

}
#endif
