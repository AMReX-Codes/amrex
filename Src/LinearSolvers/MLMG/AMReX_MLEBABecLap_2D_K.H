#ifndef AMREX_MLEBABECLAP_2D_K_H_
#define AMREX_MLEBABECLAP_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_adotx (Box const& box, Array4<Real> const& y,
                        Array4<Real const> const& x, Array4<Real const> const& a,
                        Array4<Real const> const& bX, Array4<Real const> const& bY,
                        Array4<const int> const& ccm, Array4<EBCellFlag const> const& flag,
                        Array4<Real const> const& vfrc, Array4<Real const> const& apx,
                        Array4<Real const> const& apy, Array4<Real const> const& fcx,
                        Array4<Real const> const& fcy, Array4<Real const> const& ba,
                        Array4<Real const> const& bc, Array4<Real const> const& beb,
                        bool is_dirichlet, Array4<Real const> const& phieb,
                        bool is_inhomog, GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                        Real alpha, Real beta, int ncomp)
{
    Real dhx = beta*dxinv[0]*dxinv[0];
    Real dhy = beta*dxinv[1]*dxinv[1];

    amrex::Loop(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isCovered())
        {
            y(i,j,k,n) = 0.0;
        }
        else if (flag(i,j,k).isRegular())
        {
            y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n)
                - dhx * (bX(i+1,j,k,n)*(x(i+1,j,k,n) - x(i  ,j,k,n))
                       - bX(i  ,j,k,n)*(x(i  ,j,k,n) - x(i-1,j,k,n)))
                - dhy * (bY(i,j+1,k,n)*(x(i,j+1,k,n) - x(i,j  ,k,n))
                       - bY(i,j  ,k,n)*(x(i,j  ,k,n) - x(i,j-1,k,n)));
        }
        else
        {
            Real kappa = vfrc(i,j,k);
            Real apxm = apx(i,j,k);
            Real apxp = apx(i+1,j,k);
            Real apym = apy(i,j,k);
            Real apyp = apy(i,j+1,k);

            Real fxm = bX(i,j,k,n)*(x(i,j,k,n)-x(i-1,j,k,n));
            if (apxm != 0.0 and apxm != 1.0) {
                int jj = j + static_cast<int>(std::copysign(1.0,fcx(i,j,k)));
                Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx(i,j,k)) : 0.0;
                fxm = (1.0-fracy)*fxm + fracy*bX(i,jj,k,n)*(x(i,jj,k,n)-x(i-1,jj,k,n));
            }

            Real fxp = bX(i+1,j,k,n)*(x(i+1,j,k,n)-x(i,j,k,n));
            if (apxp != 0.0 and apxp != 1.0) {
                int jj = j + static_cast<int>(std::copysign(1.0,fcx(i+1,j,k)));
                Real fracy = (ccm(i,jj,k) || ccm(i+1,jj,k)) ? std::abs(fcx(i+1,j,k)) : 0.0;
                fxp = (1.0-fracy)*fxp + fracy*bX(i+1,jj,k,n)*(x(i+1,jj,k,n)-x(i,jj,k,n));
            }

            Real fym = bY(i,j,k,n)*(x(i,j,k,n)-x(i,j-1,k,n));
            if (apym != 0.0 and apym != 1.0) {
                int ii = i + static_cast<int>(std::copysign(1.0,fcy(i,j,k)));
                Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy(i,j,k)) : 0.0;
                fym = (1.0-fracx)*fym + fracx*bY(ii,j,k,n)*(x(ii,j,k,n)-x(ii,j-1,k,n));
            }

            Real fyp = bY(i,j+1,k,n)*(x(i,j+1,k,n)-x(i,j,k,n));
            if (apyp != 0.0 and apyp != 1.0) {
                int ii = i + static_cast<int>(std::copysign(1.0,fcy(i,j+1,k)));
                Real fracx = (ccm(ii,j,k) || ccm(ii,j+1,k)) ? std::abs(fcy(i,j+1,k)) : 0.0;
                fyp = (1.0-fracx)*fyp + fracx*bY(ii,j+1,k,n)*(x(ii,j+1,k,n)-x(ii,j,k,n));
            }

            Real feb = 0.0;
            if (is_dirichlet) {
                Real dapx = apxm-apxp;
                Real dapy = apym-apyp;
                Real anorm = std::hypot(dapx,dapy);
                Real anorminv = 1.0/anorm;
                Real anrmx = dapx * anorminv;
                Real anrmy = dapy * anorminv;

                Real phib = is_inhomog ? phieb(i,j,k,n) : 0.0;

                Real bctx = bc(i,j,k,0);
                Real bcty = bc(i,j,k,1);
                Real dx_eb = get_dx_eb(kappa);

                Real dg, gx, gy, sx, sy;
                if (std::abs(anrmx) > std::abs(anrmy)) {
                    dg = dx_eb / std::abs(anrmx);
                    gx = bctx - dg*anrmx;
                    gy = bcty - dg*anrmy;
                    sx = std::copysign(1.0,anrmx);
                    sy = std::copysign(1.0,anrmy);
                } else {
                    dg = dx_eb / std::abs(anrmy);
                    gx = bctx - dg*anrmx;
                    gy = bcty - dg*anrmy;
                    sx = std::copysign(1.0,anrmx);
                    sy = std::copysign(1.0,anrmy);
                }

                int ii = i - static_cast<int>(sx);
                int jj = j - static_cast<int>(sy);

                Real phig = (1.0 + gx*sx + gy*sy + gx*gy*sx*sy) * x(i ,j ,k,n)
                    +       (    - gx*sx         - gx*gy*sx*sy) * x(ii,j ,k,n)
                    +       (            - gy*sy - gx*gy*sx*sy) * x(i ,jj,k,n)
                    +       (                    + gx*gy*sx*sy) * x(ii,jj,k,n) ;

                Real dphidn = (phib-phig) / dg;

                feb = dphidn * ba(i,j,k) * beb(i,j,k,n);
            }

            y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n) + (1.0/kappa) *
                (dhx*(apxm*fxm-apxp*fxp) +
                 dhy*(apym*fym-apyp*fyp) - dhx*feb);
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_gsrb (Box const& box,
                       Array4<Real> const& phi, Array4<Real const> const& rhs,
                       Real alpha, Array4<Real const> const& a,
                       Real dhx, Real dhy,
                       Array4<Real const> const& bX, Array4<Real const> const& bY,
                       Array4<int const> const& m0, Array4<int const> const& m2,
                       Array4<int const> const& m1, Array4<int const> const& m3,
                       Array4<Real const> const& f0, Array4<Real const> const& f2,
                       Array4<Real const> const& f1, Array4<Real const> const& f3, 
                       Array4<const int> const& ccm, Array4<EBCellFlag const> const& flag,
                       Array4<Real const> const& vfrc,
                       Array4<Real const> const& apx, Array4<Real const> const& apy,
                       Array4<Real const> const& fcx, Array4<Real const> const& fcy,
                       Array4<Real const> const& ba, Array4<Real const> const& bc,
                       Array4<Real const> const& beb,
                       bool is_dirichlet, Box const& vbox, int redblack, int ncomp)
{
    const auto vlo = amrex::lbound(vbox);
    const auto vhi = amrex::ubound(vbox);

    amrex::Loop(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if ((i+j+k+redblack) % 2 == 0)
        {
            if (flag(i,j,k).isCovered())
            {
                phi(i,j,k,n) = 0.0;
            }
            else
            {
                Real cf0 = (i == vlo.x and m0(vlo.x-1,j,k) > 0)
                    ? f0(vlo.x,j,k,n) : 0.0;
                Real cf1 = (j == vlo.y and m1(i,vlo.y-1,k) > 0)
                    ? f1(i,vlo.y,k,n) : 0.0;
                Real cf2 = (i == vhi.x and m2(vhi.x+1,j,k) > 0)
                    ? f2(vhi.x,j,k,n) : 0.0;
                Real cf3 = (j == vhi.y and m3(i,vhi.y+1,k) > 0)
                    ? f3(i,vhi.y,k,n) : 0.0;

                if (flag(i,j,k).isRegular())
                {
                    Real gamma = alpha*a(i,j,k)
                        + dhx * (bX(i+1,j,k,n) + bX(i,j,k,n))
                        + dhy * (bY(i,j+1,k,n) + bY(i,j,k,n));

                    Real rho =  dhx * (bX(i+1,j,k,n)*phi(i+1,j,k,n)
                                     + bX(i  ,j,k,n)*phi(i-1,j,k,n))
                              + dhy * (bY(i,j+1,k,n)*phi(i,j+1,k,n)
                                     + bY(i,j  ,k,n)*phi(i,j-1,k,n));

                    Real delta = dhx*(bX(i,j,k,n)*cf0 + bX(i+1,j,k,n)*cf2)
                        +        dhy*(bY(i,j,k,n)*cf1 + bY(i,j+1,k,n)*cf3);

                    Real res = rhs(i,j,k,n) - (gamma*phi(i,j,k,n) - rho);
                    phi(i,j,k,n) += res/(gamma-delta);
                }
                else
                {
                    Real kappa = vfrc(i,j,k);
                    Real apxm = apx(i,j,k);
                    Real apxp = apx(i+1,j,k);
                    Real apym = apy(i,j,k);
                    Real apyp = apy(i,j+1,k);

                    Real fxm = -bX(i,j,k,n)*phi(i-1,j,k,n);
                    Real oxm = -bX(i,j,k,n)*cf0;
                    Real sxm =  bX(i,j,k,n);
                    if (apxm != 0.0 and apxm != 1.0) {
                        int jj = j + static_cast<int>(std::copysign(1.0,fcx(i,j,k)));
                        Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k))
                            ? std::abs(fcx(i,j,k)) : 0.0;
                        fxm = (1.0-fracy)*fxm + fracy*bX(i,jj,k,n)*
                            (phi(i,jj,k,n)-phi(i-1,jj,k,n));
                        oxm = 0.0;
                        sxm = (1.0-fracy)*sxm;
                    }

                    Real fxp =  bX(i+1,j,k,n)*phi(i+1,j,k,n);
                    Real oxp =  bX(i+1,j,k,n)*cf2;
                    Real sxp = -bX(i+1,j,k,n);
                    if (apxp != 0.0 and apxp != 1.0) {
                        int jj = j + static_cast<int>(std::copysign(1.0,fcx(i+1,j,k)));
                        Real fracy = (ccm(i,jj,k) || ccm(i+1,jj,k))
                            ? std::abs(fcx(i+1,j,k)) : 0.0;
                        fxp = (1.0-fracy)*fxp + fracy*bX(i+1,jj,k,n)*
                            (phi(i+1,jj,k,n)-phi(i,jj,k,n));
                        oxp = 0.0;
                        sxp = (1.0-fracy)*sxp;
                    }

                    Real fym = -bY(i,j,k,n)*phi(i,j-1,k,n);
                    Real oym = -bY(i,j,k,n)*cf1;
                    Real sym =  bY(i,j,k,n);
                    if (apym != 0.0 and apym != 1.0) {
                        int ii = i + static_cast<int>(std::copysign(1.0,fcy(i,j,k)));
                        Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k))
                            ? std::abs(fcy(i,j,k)) : 0.0;
                        fym = (1.0-fracx)*fym + fracx*bY(ii,j,k,n)*
                            (phi(ii,j,k,n)-phi(ii,j-1,k,n));
                        oym = 0.0;
                        sym = (1.0-fracx)*sym;
                    }

                    Real fyp =  bY(i,j+1,k,n)*phi(i,j+1,k,n);
                    Real oyp =  bY(i,j+1,k,n)*cf3;
                    Real syp = -bY(i,j+1,k,n);
                    if (apyp != 0.0 and apyp != 1.0) {
                        int ii = i + static_cast<int>(std::copysign(1.0,fcy(i,j+1,k)));
                        Real fracx = (ccm(ii,j,k) || ccm(ii,j+1,k))
                            ? std::abs(fcy(i,j+1,k)) : 0.0;
                        fyp = (1.0-fracx)*fyp + fracx*bY(ii,j+1,k,n)*
                            (phi(ii,j+1,k,n)-phi(ii,j,k,n));
                        oyp = 0.0;
                        syp = (1.0-fracx)*syp;
                    }

                    Real vfrcinv = (1.0/kappa);
                    Real gamma = alpha*a(i,j,k) + vfrcinv *
                        (dhx*(apxm*sxm-apxp*sxp) +
                         dhy*(apym*sym-apyp*syp));
                    Real rho = -vfrcinv *
                        (dhx*(apxm*fxm-apxp*fxp) +
                         dhy*(apym*fym-apyp*fyp));

                    Real delta = -vfrcinv *
                        (dhx*(apxm*oxm-apxp*oxp) +
                         dhy*(apym*oym-apyp*oyp));

                    if (is_dirichlet) {
                        Real dapx = apxm-apxp;
                        Real dapy = apym-apyp;
                        Real anorm = std::hypot(dapx,dapy);
                        Real anorminv = 1.0/anorm;
                        Real anrmx = dapx * anorminv;
                        Real anrmy = dapy * anorminv;

                        Real bctx = bc(i,j,k,0);
                        Real bcty = bc(i,j,k,1);
                        Real dx_eb = get_dx_eb(vfrc(i,j,k));

                        Real dg, gx, gy, sx, sy;
                        if (std::abs(anrmx) > std::abs(anrmy)) {
                            dg = dx_eb / std::abs(anrmx);
                            gx = bctx - dg*anrmx;
                            gy = bcty - dg*anrmy;
                            sx = std::copysign(1.0,anrmx);
                            sy = std::copysign(1.0,anrmy);
                        } else {
                            dg = dx_eb / std::abs(anrmy);
                            gx = bctx - dg*anrmx;
                            gy = bcty - dg*anrmy;
                            sx = std::copysign(1.0,anrmx);
                            sy = std::copysign(1.0,anrmy);
                        }

                        int ii = i - static_cast<int>(sx);
                        int jj = j - static_cast<int>(sy);
                   
                        Real phig_gamma = (1.0 + gx*sx + gy*sy + gx*gy*sx*sy);
                        Real phig = (    - gx*sx         - gx*gy*sx*sy) * phi(ii,j ,k,n)
                            +       (            - gy*sy - gx*gy*sx*sy) * phi(i ,jj,k,n)
                            +       (                    + gx*gy*sx*sy) * phi(ii,jj,k,n);

                        // In gsrb we are always in residual-correction form so phib = 0
                        Real dphidn =  (    -phig)/dg;

                        Real feb = dphidn * ba(i,j,k) * beb(i,j,k,n);
                        rho += -vfrcinv*(-dhx)*feb;

                        Real feb_gamma = -phig_gamma/dg * ba(i,j,k) * beb(i,j,k,n);
                        gamma += vfrcinv*(-dhx)*feb_gamma;
                    }

                    Real res = rhs(i,j,k,n) - (gamma*phi(i,j,k,n) - rho);
                    phi(i,j,k,n) += res/(gamma-delta);
                }
            }
        }
    });
}

}

#endif
