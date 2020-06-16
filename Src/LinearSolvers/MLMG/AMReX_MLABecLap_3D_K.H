#ifndef AMREX_MLABECLAP_3D_K_H_
#define AMREX_MLABECLAP_3D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_adotx (Box const& box, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Array4<Real const> const& a,
                      Array4<Real const> const& bX,
                      Array4<Real const> const& bY,
                      Array4<Real const> const& bZ,
                      GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                      Real alpha, Real beta, int ncomp) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];
    const Real dhz = beta*dxinv[2]*dxinv[2];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n)
                    - dhx * (bX(i+1,j,k,n)*(x(i+1,j,k,n) - x(i  ,j,k,n))
                           - bX(i  ,j,k,n)*(x(i  ,j,k,n) - x(i-1,j,k,n)))
                    - dhy * (bY(i,j+1,k,n)*(x(i,j+1,k,n) - x(i,j  ,k,n))
                           - bY(i,j  ,k,n)*(x(i,j  ,k,n) - x(i,j-1,k,n)))
                    - dhz * (bZ(i,j,k+1,n)*(x(i,j,k+1,n) - x(i,j,k  ,n))
                           - bZ(i,j,k  ,n)*(x(i,j,k  ,n) - x(i,j,k-1,n)));
            }
        }
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_adotx_os (Box const& box, Array4<Real> const& y,
                         Array4<Real const> const& x,
                         Array4<Real const> const& a,
                         Array4<Real const> const& bX,
                         Array4<Real const> const& bY,
                         Array4<Real const> const& bZ,
                         Array4<int const> const& osm,
                         GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                         Real alpha, Real beta, int ncomp) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];
    const Real dhz = beta*dxinv[2]*dxinv[2];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (osm(i,j,k)) {
                    y(i,j,k,n) = 0.0;
                } else {
                    y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n)
                        - dhx * (bX(i+1,j,k,n)*(x(i+1,j,k,n) - x(i  ,j,k,n))
                               - bX(i  ,j,k,n)*(x(i  ,j,k,n) - x(i-1,j,k,n)))
                        - dhy * (bY(i,j+1,k,n)*(x(i,j+1,k,n) - x(i,j  ,k,n))
                               - bY(i,j  ,k,n)*(x(i,j  ,k,n) - x(i,j-1,k,n)))
                        - dhz * (bZ(i,j,k+1,n)*(x(i,j,k+1,n) - x(i,j,k  ,n))
                               - bZ(i,j,k  ,n)*(x(i,j,k  ,n) - x(i,j,k-1,n)));
                }
            }
        }
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_normalize (Box const& box, Array4<Real> const& x,
                          Array4<Real const> const& a,
                          Array4<Real const> const& bX,
                          Array4<Real const> const& bY,
                          Array4<Real const> const& bZ,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                          Real alpha, Real beta, int ncomp) noexcept
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];
    const Real dhz = beta*dxinv[2]*dxinv[2];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                x(i,j,k,n) /= alpha*a(i,j,k)
                    + dhx*(bX(i,j,k,n)+bX(i+1,j,k,n))
                    + dhy*(bY(i,j,k,n)+bY(i,j+1,k,n))
                    + dhz*(bZ(i,j,k,n)+bZ(i,j,k+1,n));
            }
        }
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_x (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                       Array4<Real const> const& bx, Real fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fx(i,j,k,n) = -fac*bx(i,j,k,n)*(sol(i,j,k,n)-sol(i-1,j,k,n));
            }
        }
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_xface (Box const& box, Array4<Real> const& fx, Array4<Real const> const& sol,
                           Array4<Real const> const& bx, Real fac, int xlen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            int i = lo.x;
            fx(i,j,k,n) = -fac*bx(i,j,k,n)*(sol(i,j,k,n)-sol(i-1,j,k,n));
            i += xlen;
            fx(i,j,k,n) = -fac*bx(i,j,k,n)*(sol(i,j,k,n)-sol(i-1,j,k,n));
        }
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_y (Box const& box, Array4<Real> const& fy, Array4<Real const> const& sol,
                       Array4<Real const> const& by, Real fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fy(i,j,k,n) = -fac*by(i,j,k,n)*(sol(i,j,k,n)-sol(i,j-1,k,n));
            }
        }
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_yface (Box const& box, Array4<Real> const& fy, Array4<Real const> const& sol,
                           Array4<Real const> const& by, Real fac, int ylen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for     (int k = lo.z; k <= hi.z; ++k) {
        int j = lo.y;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,k,n) = -fac*by(i,j,k,n)*(sol(i,j,k,n)-sol(i,j-1,k,n));
        }
        j += ylen;
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,k,n) = -fac*by(i,j,k,n)*(sol(i,j,k,n)-sol(i,j-1,k,n));
        }
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_z (Box const& box, Array4<Real> const& fz, Array4<Real const> const& sol,
                       Array4<Real const> const& bz, Real fac, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                fz(i,j,k,n) = -fac*bz(i,j,k,n)*(sol(i,j,k,n)-sol(i,j,k-1,n));
            }
        }
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlabeclap_flux_zface (Box const& box, Array4<Real> const& fz, Array4<Real const> const& sol,
                           Array4<Real const> const& bz, Real fac, int zlen, int ncomp) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < ncomp; ++n) {
    int k = lo.z;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fz(i,j,k,n) = -fac*bz(i,j,k,n)*(sol(i,j,k,n)-sol(i,j,k-1,n));
        }
    }

    k += zlen;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fz(i,j,k,n) = -fac*bz(i,j,k,n)*(sol(i,j,k,n)-sol(i,j,k-1,n));
        }
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void abec_gsrb (Box const& box, Array4<Real> const& phi, Array4<Real const> const& rhs,
                Real alpha, Array4<Real const> const& a,
                Real dhx, Real dhy, Real dhz,
                Array4<Real const> const& bX, Array4<Real const> const& bY,
                Array4<Real const> const& bZ,
                Array4<int const> const& m0, Array4<int const> const& m2,
                Array4<int const> const& m4,
                Array4<int const> const& m1, Array4<int const> const& m3,
                Array4<int const> const& m5,
                Array4<Real const> const& f0, Array4<Real const> const& f2,
                Array4<Real const> const& f4,
                Array4<Real const> const& f1, Array4<Real const> const& f3,
                Array4<Real const> const& f5,
                Box const& vbox, int redblack, int nc) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    constexpr Real omega = 1.15;

    for (int n = 0; n < nc; ++n) {
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    if ((i+j+k+redblack)%2 == 0) {
                        Real cf0 = (i == vlo.x and m0(vlo.x-1,j,k) > 0)
                            ? f0(vlo.x,j,k,n) : 0.0;
                        Real cf1 = (j == vlo.y and m1(i,vlo.y-1,k) > 0)
                            ? f1(i,vlo.y,k,n) : 0.0;
                        Real cf2 = (k == vlo.z and m2(i,j,vlo.z-1) > 0)
                            ? f2(i,j,vlo.z,n) : 0.0;
                        Real cf3 = (i == vhi.x and m3(vhi.x+1,j,k) > 0)
                            ? f3(vhi.x,j,k,n) : 0.0;
                        Real cf4 = (j == vhi.y and m4(i,vhi.y+1,k) > 0)
                            ? f4(i,vhi.y,k,n) : 0.0;
                        Real cf5 = (k == vhi.z and m5(i,j,vhi.z+1) > 0)
                            ? f5(i,j,vhi.z,n) : 0.0;

                        Real gamma = alpha*a(i,j,k)
                            +   dhx*(bX(i,j,k,n)+bX(i+1,j,k,n))
                            +   dhy*(bY(i,j,k,n)+bY(i,j+1,k,n))
                            +   dhz*(bZ(i,j,k,n)+bZ(i,j,k+1,n));

                        Real g_m_d = gamma
                            - (dhx*(bX(i,j,k,n)*cf0 + bX(i+1,j,k,n)*cf3)
                            +  dhy*(bY(i,j,k,n)*cf1 + bY(i,j+1,k,n)*cf4)
                            +  dhz*(bZ(i,j,k,n)*cf2 + bZ(i,j,k+1,n)*cf5));

                        Real rho =  dhx*( bX(i  ,j,k,n)*phi(i-1,j,k,n)
                                  +       bX(i+1,j,k,n)*phi(i+1,j,k,n) )
                                  + dhy*( bY(i,j  ,k,n)*phi(i,j-1,k,n)
                                  +       bY(i,j+1,k,n)*phi(i,j+1,k,n) )
                                  + dhz*( bZ(i,j,k  ,n)*phi(i,j,k-1,n)
                                  +       bZ(i,j,k+1,n)*phi(i,j,k+1,n) );

                        Real res =  rhs(i,j,k,n) - (gamma*phi(i,j,k,n) - rho);
                        phi(i,j,k,n) = phi(i,j,k,n) + omega/g_m_d * res;
                    }
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void abec_gsrb_os (Box const& box, Array4<Real> const& phi, Array4<Real const> const& rhs,
                   Real alpha, Array4<Real const> const& a,
                   Real dhx, Real dhy, Real dhz,
                   Array4<Real const> const& bX, Array4<Real const> const& bY,
                   Array4<Real const> const& bZ,
                   Array4<int const> const& m0, Array4<int const> const& m2,
                   Array4<int const> const& m4,
                   Array4<int const> const& m1, Array4<int const> const& m3,
                   Array4<int const> const& m5,
                   Array4<Real const> const& f0, Array4<Real const> const& f2,
                   Array4<Real const> const& f4,
                   Array4<Real const> const& f1, Array4<Real const> const& f3,
                   Array4<Real const> const& f5,
                   Array4<int const> const& osm,
                   Box const& vbox, int redblack, int nc) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    constexpr Real omega = 1.15;

    for (int n = 0; n < nc; ++n) {
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    if ((i+j+k+redblack)%2 == 0) {
                        if (osm(i,j,k)) {
                            phi(i,j,k,n) = 0.0;
                        } else {
                            Real cf0 = (i == vlo.x and m0(vlo.x-1,j,k) > 0)
                                ? f0(vlo.x,j,k,n) : 0.0;
                            Real cf1 = (j == vlo.y and m1(i,vlo.y-1,k) > 0)
                                ? f1(i,vlo.y,k,n) : 0.0;
                            Real cf2 = (k == vlo.z and m2(i,j,vlo.z-1) > 0)
                                ? f2(i,j,vlo.z,n) : 0.0;
                            Real cf3 = (i == vhi.x and m3(vhi.x+1,j,k) > 0)
                                ? f3(vhi.x,j,k,n) : 0.0;
                            Real cf4 = (j == vhi.y and m4(i,vhi.y+1,k) > 0)
                                ? f4(i,vhi.y,k,n) : 0.0;
                            Real cf5 = (k == vhi.z and m5(i,j,vhi.z+1) > 0)
                                ? f5(i,j,vhi.z,n) : 0.0;

                            Real gamma = alpha*a(i,j,k)
                                +   dhx*(bX(i,j,k,n)+bX(i+1,j,k,n))
                                +   dhy*(bY(i,j,k,n)+bY(i,j+1,k,n))
                                +   dhz*(bZ(i,j,k,n)+bZ(i,j,k+1,n));

                            Real g_m_d = gamma
                                - (dhx*(bX(i,j,k,n)*cf0 + bX(i+1,j,k,n)*cf3)
                                +  dhy*(bY(i,j,k,n)*cf1 + bY(i,j+1,k,n)*cf4)
                                +  dhz*(bZ(i,j,k,n)*cf2 + bZ(i,j,k+1,n)*cf5));

                            Real rho =  dhx*( bX(i  ,j,k,n)*phi(i-1,j,k,n)
                                      +       bX(i+1,j,k,n)*phi(i+1,j,k,n) )
                                      + dhy*( bY(i,j  ,k,n)*phi(i,j-1,k,n)
                                      +       bY(i,j+1,k,n)*phi(i,j+1,k,n) )
                                      + dhz*( bZ(i,j,k  ,n)*phi(i,j,k-1,n)
                                      +       bZ(i,j,k+1,n)*phi(i,j,k+1,n) );

                            Real res =  rhs(i,j,k,n) - (gamma*phi(i,j,k,n) - rho);
                            phi(i,j,k,n) = phi(i,j,k,n) + omega/g_m_d * res;
                        }
                    }
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int coarsen_overset_mask (Box const& bx, Array4<int> const& cmsk, Array4<int const> const& fmsk)
    noexcept
{
    int nerrors = 0;
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        int ii = 2*i;
        int jj = 2*j;
        int kk = 2*k;
        cmsk(i,j,k) = fmsk(ii,jj  ,kk  ) + fmsk(ii+1,jj  ,kk  )
            +         fmsk(ii,jj+1,kk  ) + fmsk(ii+1,jj+1,kk  )
            +         fmsk(ii,jj  ,kk+1) + fmsk(ii+1,jj  ,kk+1)
            +         fmsk(ii,jj+1,kk+1) + fmsk(ii+1,jj+1,kk+1);
        if (cmsk(i,j,k) == 8) {
            cmsk(i,j,k) = 1;
        } else if (cmsk(i,j,k) != 0) {
            ++nerrors;
        }
    }}}
    return nerrors;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void overset_rescale_bcoef_x (Box const& box, Array4<Real> const& bX, Array4<int const> const& osm,
                              int ncomp, Real osfac) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((osm(i-1,j,k)+osm(i,j,k)) == 1) {
                bX(i,j,k,n) *= osfac;
            }
        }}}
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void overset_rescale_bcoef_y (Box const& box, Array4<Real> const& bY, Array4<int const> const& osm,
                              int ncomp, Real osfac) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((osm(i,j-1,k)+osm(i,j,k)) == 1) {
                bY(i,j,k,n) *= osfac;
            }
        }}}
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void overset_rescale_bcoef_z (Box const& box, Array4<Real> const& bZ, Array4<int const> const& osm,
                              int ncomp, Real osfac) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((osm(i,j,k-1)+osm(i,j,k)) == 1) {
                bZ(i,j,k,n) *= osfac;
            }
        }}}
    }
}

}
#endif
