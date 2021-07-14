#ifndef AMREX_MLALAP_3D_K_H_
#define AMREX_MLALAP_3D_K_H_
#include <AMReX_Config.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_adotx (Box const& box, Array4<Real> const& y,
                   Array4<Real const> const& x,
                   Array4<Real const> const& a,
                   GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                   Real alpha, Real beta, int ncomp) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];
    const Real dhz = beta*dxinv[2]*dxinv[2];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++ncomp) {
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    y(i,j,k,n) = alpha*a(i,j,k,n)*x(i,j,k,n)
                        - dhx * (x(i-1,j,k,n) - Real(2.0)*x(i,j,k,n) + x(i+1,j,k,n))
                        - dhy * (x(i,j-1,k,n) - Real(2.0)*x(i,j,k,n) + x(i,j+1,k,n))
                        - dhz * (x(i,j,k-1,n) - Real(2.0)*x(i,j,k,n) + x(i,j,k+1,n));
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_normalize (Box const& box, Array4<Real> const& x,
                       Array4<Real const> const& a,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                       Real alpha, Real beta, int ncomp) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];
    const Real dhz = beta*dxinv[2]*dxinv[2];
    const Real fac = Real(2.0)*(dhx+dhy+dhz);

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++ncomp) {
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    x(i,j,k,n) /= alpha*a(i,j,k,n) + fac;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_x (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                    Real fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++ncomp) {
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    fx(i,j,k,n) = -fac*(sol(i,j,k,n)-sol(i-1,j,k,n));
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_xface (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                        Real fac, int xlen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++ncomp) {
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                int i = lo.x;
                fx(i,j,k,n) = -fac*(sol(i,j,k,n)-sol(i-1,j,k,n));
                i += xlen;
                fx(i,j,k,n) = -fac*(sol(i,j,k,n)-sol(i-1,j,k,n));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_y (Box const& box, Array4<Real> const& fy, Array4<Real const> const& sol,
                    Real fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++ncomp) {
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    fy(i,j,k,n) = -fac*(sol(i,j,k,n)-sol(i,j-1,k,n));
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_yface (Box const& box, Array4<Real> const& fy, Array4<Real const> const& sol,
                        Real fac, int ylen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++ncomp) {
        for     (int k = lo.z; k <= hi.z; ++k) {
            int j = lo.y;
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fy(i,j,k,n) = -fac*(sol(i,j,k,n)-sol(i,j-1,k,n));
            }
            j += ylen;
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fy(i,j,k,n) = -fac*(sol(i,j,k,n)-sol(i,j-1,k,n));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_z (Box const& box, Array4<Real> const& fz, Array4<Real const> const& sol,
                    Real fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++ncomp) {
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    fz(i,j,k,n) = -fac*(sol(i,j,k,n)-sol(i,j,k-1,n));
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_flux_zface (Box const& box, Array4<Real> const& fz, Array4<Real const> const& sol,
                        Real fac, int zlen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++ncomp) {
        int k = lo.z;
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fz(i,j,k,n) = -fac*(sol(i,j,k,n)-sol(i,j,k-1,n));
            }
        }

        k += zlen;
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fz(i,j,k,n) = -fac*(sol(i,j,k,n)-sol(i,j,k-1,n));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlalap_gsrb (Box const& box, Array4<Real> const& phi,
                  Array4<Real const> const& rhs, Real alpha,
                  Real dhx, Real dhy, Real dhz, Array4<Real const> const& a,
                  Array4<Real const> const& f0, Array4<int const> const& m0,
                  Array4<Real const> const& f1, Array4<int const> const& m1,
                  Array4<Real const> const& f2, Array4<int const> const& m2,
                  Array4<Real const> const& f3, Array4<int const> const& m3,
                  Array4<Real const> const& f4, Array4<int const> const& m4,
                  Array4<Real const> const& f5, Array4<int const> const& m5,
                  Box const& vbox, int redblack, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    constexpr Real omega = Real(1.15);

    const Real dhfac = Real(2.)*(dhx+dhy+dhz);

    for (int n = 0; n < ncomp; ++ncomp) {
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    if ((i+j+k+redblack)%2 == 0) {
                        Real cf0 = (i == vlo.x && m0(vlo.x-1,j,k) > 0)
                            ? f0(vlo.x,j,k,n) : Real(0.0);
                        Real cf1 = (j == vlo.y && m1(i,vlo.y-1,k) > 0)
                            ? f1(i,vlo.y,k,n) : Real(0.0);
                        Real cf2 = (k == vlo.z && m2(i,j,vlo.z-1) > 0)
                            ? f2(i,j,vlo.z,n) : Real(0.0);
                        Real cf3 = (i == vhi.x && m3(vhi.x+1,j,k) > 0)
                            ? f3(vhi.x,j,k,n) : Real(0.0);
                        Real cf4 = (j == vhi.y && m4(i,vhi.y+1,k) > 0)
                            ? f4(i,vhi.y,k,n) : Real(0.0);
                        Real cf5 = (k == vhi.z && m5(i,j,vhi.z+1) > 0)
                            ? f5(i,j,vhi.z,n) : Real(0.0);

                        Real gamma = alpha*a(i,j,k,n) + dhfac;

                        Real g_m_d = gamma - dhx*(cf0+cf3) - dhy*(cf1+cf4) - dhz*(cf2+cf5);

                        Real rho =  dhx*(phi(i-1,j,k,n) + phi(i+1,j,k,n))
                                  + dhy*(phi(i,j-1,k,n) + phi(i,j+1,k,n))
                                  + dhz*(phi(i,j,k-1,n) + phi(i,j,k+1,n));

                        Real res =  rhs(i,j,k,n) - (gamma*phi(i,j,k,n) - rho);
                        phi(i,j,k,n) = phi(i,j,k,n) + omega/g_m_d * res;
                    }
                }
            }
        }
    }
}

}
#endif
