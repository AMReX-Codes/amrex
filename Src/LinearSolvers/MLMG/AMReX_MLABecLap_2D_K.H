#ifndef AMREX_MLABECLAP_2D_K_H_
#define AMREX_MLABECLAP_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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
                - dhx * (bX(i+1,j,0,n)*(x(i+1,j,0,n) - x(i  ,j,0,n))
                       - bX(i  ,j,0,n)*(x(i  ,j,0,n) - x(i-1,j,0,n)))
                - dhy * (bY(i,j+1,0,n)*(x(i,j+1,0,n) - x(i,j  ,0,n))
                       - bY(i,j  ,0,n)*(x(i,j  ,0,n) - x(i,j-1,0,n)));
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
                         Array4<int const> const& osm,
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
            if (osm(i,j,0) == 0) {
                y(i,j,0,n) = 0.0;
            } else {
                y(i,j,0,n) = alpha*a(i,j,0)*x(i,j,0,n)
                    - dhx * (bX(i+1,j,0,n)*(x(i+1,j,0,n) - x(i  ,j,0,n))
                           - bX(i  ,j,0,n)*(x(i  ,j,0,n) - x(i-1,j,0,n)))
                    - dhy * (bY(i,j+1,0,n)*(x(i,j+1,0,n) - x(i,j  ,0,n))
                           - bY(i,j  ,0,n)*(x(i,j  ,0,n) - x(i,j-1,0,n)));
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
                + dhx*(bX(i,j,0,n)+bX(i+1,j,0,n))
                + dhy*(bY(i,j,0,n)+bY(i,j+1,0,n));
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
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fx(i,j,0,n) = -fac*bx(i,j,0,n)*(sol(i,j,0,n)-sol(i-1,j,0,n));
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
    for     (int j = lo.y; j <= hi.y; ++j) {
        int i = lo.x;
        fx(i,j,0,n) = -fac*bx(i,j,0,n)*(sol(i,j,0,n)-sol(i-1,j,0,n));
        i += xlen;
        fx(i,j,0,n) = -fac*bx(i,j,0,n)*(sol(i,j,0,n)-sol(i-1,j,0,n));
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
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            fy(i,j,0,n) = -fac*by(i,j,0,n)*(sol(i,j,0,n)-sol(i,j-1,0,n));
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
    int j = lo.y;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fy(i,j,0,n) = -fac*by(i,j,0,n)*(sol(i,j,0,n)-sol(i,j-1,0,n));
    }
    j += ylen;
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        fy(i,j,0,n) = -fac*by(i,j,0,n)*(sol(i,j,0,n)-sol(i,j-1,0,n));
    }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void abec_gsrb (Box const& box, Array4<Real> const& phi, Array4<Real const> const& rhs,
                Real alpha, Array4<Real const> const& a,
                Real dhx, Real dhy,
                Array4<Real const> const& bX, Array4<Real const> const& bY,
                Array4<int const> const& m0, Array4<int const> const& m2,
                Array4<int const> const& m1, Array4<int const> const& m3,
                Array4<Real const> const& f0, Array4<Real const> const& f2,
                Array4<Real const> const& f1, Array4<Real const> const& f3,
                Box const& vbox, int redblack, int nc) noexcept
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

                    Real delta = dhx*(bX(i,j,0,n)*cf0 + bX(i+1,j,0,n)*cf2)
                              +  dhy*(bY(i,j,0,n)*cf1 + bY(i,j+1,0,n)*cf3);

                    Real gamma = alpha*a(i,j,0)
                        +   dhx*( bX(i,j,0,n) + bX(i+1,j,0,n) )
                        +   dhy*( bY(i,j,0,n) + bY(i,j+1,0,n) );

                    Real rho = dhx*(bX(i  ,j  ,0,n)*phi(i-1,j  ,0,n)
                                  + bX(i+1,j  ,0,n)*phi(i+1,j  ,0,n))
                              +dhy*(bY(i  ,j  ,0,n)*phi(i  ,j-1,0,n)
                                  + bY(i  ,j+1,0,n)*phi(i  ,j+1,0,n));

                    phi(i,j,0,n) = (rhs(i,j,0,n) + rho - phi(i,j,0,n)*delta)
                        / (gamma - delta);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void abec_gsrb_os (Box const& box, Array4<Real> const& phi, Array4<Real const> const& rhs,
                   Real alpha, Array4<Real const> const& a,
                   Real dhx, Real dhy,
                   Array4<Real const> const& bX, Array4<Real const> const& bY,
                   Array4<int const> const& m0, Array4<int const> const& m2,
                   Array4<int const> const& m1, Array4<int const> const& m3,
                   Array4<Real const> const& f0, Array4<Real const> const& f2,
                   Array4<Real const> const& f1, Array4<Real const> const& f3,
                   Array4<int const> const& osm,
                   Box const& vbox, int redblack, int nc) noexcept
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
                    if (osm(i,j,0) == 0) {
                        phi(i,j,0,n) = 0.0;
                    } else {
                        Real cf0 = (i == vlo.x and m0(vlo.x-1,j,0) > 0)
                            ? f0(vlo.x,j,0,n) : 0.0;
                        Real cf1 = (j == vlo.y and m1(i,vlo.y-1,0) > 0)
                            ? f1(i,vlo.y,0,n) : 0.0;
                        Real cf2 = (i == vhi.x and m2(vhi.x+1,j,0) > 0)
                            ? f2(vhi.x,j,0,n) : 0.0;
                        Real cf3 = (j == vhi.y and m3(i,vhi.y+1,0) > 0)
                            ? f3(i,vhi.y,0,n) : 0.0;

                        Real delta = dhx*(bX(i,j,0,n)*cf0 + bX(i+1,j,0,n)*cf2)
                                  +  dhy*(bY(i,j,0,n)*cf1 + bY(i,j+1,0,n)*cf3);

                        Real gamma = alpha*a(i,j,0)
                            +   dhx*( bX(i,j,0,n) + bX(i+1,j,0,n) )
                            +   dhy*( bY(i,j,0,n) + bY(i,j+1,0,n) );

                        Real rho = dhx*(bX(i  ,j  ,0,n)*phi(i-1,j  ,0,n)
                                      + bX(i+1,j  ,0,n)*phi(i+1,j  ,0,n))
                                  +dhy*(bY(i  ,j  ,0,n)*phi(i  ,j-1,0,n)
                                      + bY(i  ,j+1,0,n)*phi(i  ,j+1,0,n));

                        phi(i,j,0,n) = (rhs(i,j,0,n) + rho - phi(i,j,0,n)*delta)
                            / (gamma - delta);
                    }
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void abec_gsrb_with_line_solve (
                Box const& box, Array4<Real> const& phi, Array4<Real const> const& rhs,
                Real alpha, Array4<Real const> const& a,
                Real dhx, Real dhy,
                Array4<Real const> const& bX, Array4<Real const> const& bY,
                Array4<int const> const& m0, Array4<int const> const& m2,
                Array4<int const> const& m1, Array4<int const> const& m3,
                Array4<Real const> const& f0, Array4<Real const> const& f2,
                Array4<Real const> const& f1, Array4<Real const> const& f3,
                Box const& vbox, int redblack, int nc) noexcept
{

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    // amrex::Print() << "GSRB LS " << box << " " << dhx << " " << dhy << " " << dhz << std::endl;

    // idir is the direction in which we will do the tridiagonal solve --
    // it should be the direction in which the mesh spacing is much larger
    // than in the other directions
    // int idir = 1;

    // This should be moved outside the kernel!
    if (dhy <= dhx) amrex::Abort("dhy is supposed to be much larger than dhx");

    int ilen = hi.y - lo.y + 1;

    // This should be moved outside the kernel!
    if (ilen > 32) amrex::Abort("abec_gsrb_with_line_solve is hard-wired to be no longer than 32");

    Array1D<Real,0,31> a_ls;
    Array1D<Real,0,31> b_ls;
    Array1D<Real,0,31> c_ls;
    Array1D<Real,0,31> r_ls;
    Array1D<Real,0,31> u_ls;
    Array1D<Real,0,31> gam;

    for (int n = 0; n < nc; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((i+redblack)%2 == 0) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    Real gamma = alpha*a(i,j,0)
                        +   dhx*(bX(i,j,0,n)+bX(i+1,j,0,n))
                        +   dhy*(bY(i,j,0,n)+bY(i,j+1,0,n));

                    Real cf0 = (i == vlo.x and m0(vlo.x-1,j,0) > 0)
                        ? f0(vlo.x,j,0,n) : 0.0;
                    Real cf1 = (j == vlo.y and m1(i,vlo.y-1,0) > 0)
                        ? f1(i,vlo.y,0,n) : 0.0;
                    Real cf2 = (i == vhi.x and m2(vhi.x+1,j,0) > 0)
                        ? f2(vhi.x,j,0,n) : 0.0;
                    Real cf3 = (j == vhi.y and m3(i,vhi.y+1,0) > 0)
                        ? f3(i,vhi.y,0,n) : 0.0;

                    Real g_m_d = gamma
                        - (dhx*(bX(i,j,0,n)*cf0 + bX(i+1,j,0,n)*cf2)
                        +  dhy*(bY(i,j,0,n)*cf1 + bY(i,j+1,0,n)*cf3));

                    Real rho =  dhx*( bX(i  ,j,0,n)*phi(i-1,j,0,n)
                              +       bX(i+1,j,0,n)*phi(i+1,j,0,n) );

                    // We have already accounted for this external boundary in the coefficient of phi(i,j,k,n)
                    if (i == vlo.x and m0(vlo.x-1,j,0) > 0)
                        rho -= dhx*bX(i  ,j,0,n)*phi(i-1,j,0,n);
                    if (i == vhi.x and m3(vhi.x+1,j,0) > 0)
                        rho -= dhx*bX(i+1,j,0,n)*phi(i+1,j,0,n);

                    a_ls(j-lo.y) = -dhy*bY(i,j,0,n);
                    b_ls(j-lo.y) =  g_m_d;
                    c_ls(j-lo.y) = -dhy*bY(i,j+1,0,n);
                    u_ls(j-lo.y) = 0.;
                    r_ls(j-lo.y) = rhs(i,j,0,n) + rho;

                    if (j == lo.y) {
                        a_ls(j-lo.y) = 0.;
                        if (!(m1(i,vlo.y-1,0) > 0)) r_ls(j-lo.y) += dhy*bY(i,j,0,n)*phi(i,j-1,0,n);
                    }
                    if (j == hi.y) {
                        c_ls(j-lo.y) = 0.;
                        if (!(m3(i,vhi.y+1,0) > 0)) r_ls(j-lo.y) += dhy*bY(i,j+1,0,n)*phi(i,j+1,0,n);
                    }
                }
//                      This is the tridiagonal solve
                {
                    Real bet = b_ls(0);
                    u_ls(0) = r_ls(0) / bet;

                    for (int jj = 1; jj <= ilen-1; jj++) {
                        gam(jj) = c_ls(jj-1) / bet;
                        bet = b_ls(jj) - a_ls(jj)*gam(jj);
                        if (bet == 0) amrex::Abort(">>>TRIDIAG FAILED");
                        u_ls(jj) = (r_ls(jj)-a_ls(jj)*u_ls(jj-1)) / bet;
                    }

                    for (int jj = ilen-2; jj >= 0; jj--) {
                                u_ls(jj) = u_ls(jj) - gam(jj+1)*u_ls(jj+1);
                    }
                }

                for (int j = lo.y; j <= hi.y; ++j) {
                            phi(i,j,0,n) = u_ls(j-lo.y);
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
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        int ii = 2*i;
        int jj = 2*j;
        cmsk(i,j,0) = fmsk(ii,jj,0) + fmsk(ii+1,jj,0) + fmsk(ii,jj+1,0) + fmsk(ii+1,jj+1,0);
        if (cmsk(i,j,0) == 4) {
            cmsk(i,j,0) = 1;
        } else if (cmsk(i,j,0) != 0) {
            ++nerrors;
        }
    }}
    return nerrors;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void overset_rescale_bcoef_x (Box const& box, Array4<Real> const& bX, Array4<int const> const& osm,
                              int ncomp, Real osfac) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int n = 0; n < ncomp; ++n) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((osm(i-1,j,0)+osm(i,j,0)) == 1) {
                bX(i,j,0,n) *= osfac;
            }
        }}
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void overset_rescale_bcoef_y (Box const& box, Array4<Real> const& bY, Array4<int const> const& osm,
                              int ncomp, Real osfac) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int n = 0; n < ncomp; ++n) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if ((osm(i,j-1,0)+osm(i,j,0)) == 1) {
                bY(i,j,0,n) *= osfac;
            }
        }}
    }
}

}
#endif
