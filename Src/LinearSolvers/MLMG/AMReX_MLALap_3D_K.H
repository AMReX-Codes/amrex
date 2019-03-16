#ifndef AMREX_MLALAP_3D_K_H_
#define AMREX_MLALAP_3D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_adotx (Box const& box, Array4<Real> const& y,
                   Array4<Real const> const& x,
                   Array4<Real const> const& a,
                   GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                   Real alpha, Real beta) noexcept
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
                    - dhx * (x(i-1,j,k) - 2.0*x(i,j,k) + x(i+1,j,k))
                    - dhy * (x(i,j-1,k) - 2.0*x(i,j,k) + x(i,j+1,k))
                    - dhz * (x(i,j,k-1) - 2.0*x(i,j,k) + x(i,j,k+1));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_normalize (Box const& box, Array4<Real> const& x,
                       Array4<Real const> const& a,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                       Real alpha, Real beta) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];
    const Real dhz = beta*dxinv[2]*dxinv[2];
    const Real fac = 2.0*(dhx+dhy+dhz);

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                x(i,j,k) /= alpha*a(i,j,k) + fac;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_flux_x (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol, Real fac) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fx(i,j,k) = -fac*(sol(i,j,k)-sol(i-1,j,k));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_flux_xface (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                        Real fac, int xlen) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            int i = lo.x;
            fx(i,j,k) = -fac*(sol(i,j,k)-sol(i-1,j,k));
            i += xlen;
            fx(i,j,k) = -fac*(sol(i,j,k)-sol(i-1,j,k));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_flux_y (Box const& box, Array4<Real> const& fy, Array4<Real const> const& sol, Real fac) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fy(i,j,k) = -fac*(sol(i,j,k)-sol(i,j-1,k));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_flux_yface (Box const& box, Array4<Real> const& fy, Array4<Real const> const& sol,
                        Real fac, int ylen) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int k = lo.z; k <= hi.z; ++k) {
        int j = lo.y;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,k) = -fac*(sol(i,j,k)-sol(i,j-1,k));
        }
        j += ylen;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,k) = -fac*(sol(i,j,k)-sol(i,j-1,k));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_flux_z (Box const& box, Array4<Real> const& fz, Array4<Real const> const& sol,
                    Real fac) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fz(i,j,k) = -fac*(sol(i,j,k)-sol(i,j,k-1));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_flux_zface (Box const& box, Array4<Real> const& fz, Array4<Real const> const& sol,
                        Real fac, int zlen) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = lo.z;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fz(i,j,k) = -fac*(sol(i,j,k)-sol(i,j,k-1));
        }
    }

    k += zlen;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fz(i,j,k) = -fac*(sol(i,j,k)-sol(i,j,k-1));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlalap_gsrb (Box const& box, Array4<Real> const& phi,
                  Array4<Real const> const& rhs, Real alpha,
                  Real dhx, Real dhy, Real dhz, Array4<Real const> const& a,
                  Array4<Real const> const& f0, Array4<int const> const& m0,
                  Array4<Real const> const& f1, Array4<int const> const& m1,
                  Array4<Real const> const& f2, Array4<int const> const& m2,
                  Array4<Real const> const& f3, Array4<int const> const& m3,
                  Array4<Real const> const& f4, Array4<int const> const& m4,
                  Array4<Real const> const& f5, Array4<int const> const& m5,
                  Box const& vbox, int redblack) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    constexpr Real omega = 1.15;

    const Real dhfac = 2.*(dhx+dhy+dhz);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if ((i+j+k+redblack)%2 == 0) {
                    Real cf0 = (i == vlo.x and m0(vlo.x-1,j,k) > 0)
                        ? f0(vlo.x,j,k) : 0.0;
                    Real cf1 = (j == vlo.y and m1(i,vlo.y-1,k) > 0)
                        ? f1(i,vlo.y,k) : 0.0;
                    Real cf2 = (k == vlo.z and m2(i,j,vlo.z-1) > 0)
                        ? f2(i,j,vlo.z) : 0.0;
                    Real cf3 = (i == vhi.x and m3(vhi.x+1,j,k) > 0)
                        ? f3(vhi.x,j,k) : 0.0;
                    Real cf4 = (j == vhi.y and m4(i,vhi.y+1,k) > 0)
                        ? f4(i,vhi.y,k) : 0.0;
                    Real cf5 = (k == vhi.z and m5(i,j,vhi.z+1) > 0)
                        ? f5(i,j,vhi.z) : 0.0;

                    Real gamma = alpha*a(i,j,k) + dhfac;

                    Real g_m_d = gamma - dhx*(cf0+cf3) - dhy*(cf1+cf4) - dhz*(cf2+cf5);

                    Real rho =  dhx*(phi(i-1,j,k) + phi(i+1,j,k))
                              + dhy*(phi(i,j-1,k) + phi(i,j+1,k))
                              + dhz*(phi(i,j,k-1) + phi(i,j,k+1));

                    Real res =  rhs(i,j,k) - (gamma*phi(i,j,k) - rho);
                    phi(i,j,k) = phi(i,j,k) + omega/g_m_d * res;
                }
            }
        }
    }
}

}
#endif
