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
                        Real alpha, Real beta, int ncomp, 
                        bool beta_on_centroid, bool phi_on_centroid) noexcept
{
    Real dhx = beta*dxinv[0]*dxinv[0];
    Real dhy = beta*dxinv[1]*dxinv[1];

    bool beta_on_center = !(beta_on_centroid);
    bool  phi_on_center = !( phi_on_centroid);

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

            Real fxm = bX(i,j,k,n) * (x(i,j,k,n)-x(i-1,j,k,n));
            if (apxm != 0.0 and apxm != 1.0) {
                int jj = j + static_cast<int>(amrex::Math::copysign(1.0,fcx(i,j,k)));
                Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? amrex::Math::abs(fcx(i,j,k)) : 0.0;
                if (beta_on_center and phi_on_center) 
                   fxm = (1.0-fracy)*fxm + fracy*bX(i,jj,k,n)*(x(i,jj,k,n)-x(i-1,jj,k,n));
                else if (beta_on_centroid and phi_on_center) 
                   fxm = bX(i,j,k,n) * ( (1.0-fracy)*(x(i, j,k,n)-x(i-1, j,k,n))
                                           + fracy  *(x(i,jj,k,n)-x(i-1,jj,k,n)) );
            }

            Real fxp = bX(i+1,j,k,n)*(x(i+1,j,k,n)-x(i,j,k,n));
            if (apxp != 0.0 and apxp != 1.0) {
                int jj = j + static_cast<int>(amrex::Math::copysign(1.0,fcx(i+1,j,k)));
                Real fracy = (ccm(i,jj,k) || ccm(i+1,jj,k)) ? amrex::Math::abs(fcx(i+1,j,k)) : 0.0;
                if (beta_on_center and phi_on_center) 
                    fxp = (1.0-fracy)*fxp + fracy*bX(i+1,jj,k,n)*(x(i+1,jj,k,n)-x(i,jj,k,n));
                else if (beta_on_centroid and phi_on_center) 
                    fxp = bX(i+1,j,k,n) * ( (1.0-fracy)*(x(i+1, j,k,n)-x(i, j,k,n))
                                               + fracy *(x(i+1,jj,k,n)-x(i,jj,k,n)) );
            }

            Real fym = bY(i,j,k,n)*(x(i,j,k,n)-x(i,j-1,k,n));
            if (apym != 0.0 and apym != 1.0) {
                int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j,k)));
                Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? amrex::Math::abs(fcy(i,j,k)) : 0.0;
                if (beta_on_center and phi_on_center) 
                    fym = (1.0-fracx)*fym + fracx*bY(ii,j,k,n)*(x(ii,j,k,n)-x(ii,j-1,k,n));
                else if (beta_on_centroid and phi_on_center) 
                    fym = bY(i,j,k,n) * ( (1.0-fracx)*(x( i,j,k,n)-x( i,j-1,k,n))
                                             + fracx *(x(ii,j,k,n)-x(ii,j-1,k,n)) );
            }

            Real fyp = bY(i,j+1,k,n)*(x(i,j+1,k,n)-x(i,j,k,n));
            if (apyp != 0.0 and apyp != 1.0) {
                int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j+1,k)));
                Real fracx = (ccm(ii,j,k) || ccm(ii,j+1,k)) ? amrex::Math::abs(fcy(i,j+1,k)) : 0.0;
                if (beta_on_center and phi_on_center) 
                    fyp = (1.0-fracx)*fyp + fracx*bY(ii,j+1,k,n)*(x(ii,j+1,k,n)-x(ii,j,k,n));
                else if (beta_on_centroid and phi_on_center) 
                    fyp = bY(i,j+1,k,n) * ( (1.0-fracx)*(x( i,j+1,k,n)-x( i,j,k,n))
                                               + fracx *(x(ii,j+1,k,n)-x(ii,j,k,n)) );
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
                if (amrex::Math::abs(anrmx) > amrex::Math::abs(anrmy)) {
                    dg = dx_eb / amrex::Math::abs(anrmx);
                    gx = bctx - dg*anrmx;
                    gy = bcty - dg*anrmy;
                    sx = amrex::Math::copysign(1.0,anrmx);
                    sy = amrex::Math::copysign(1.0,anrmy);
                } else {
                    dg = dx_eb / amrex::Math::abs(anrmy);
                    gx = bctx - dg*anrmx;
                    gy = bcty - dg*anrmy;
                    sx = amrex::Math::copysign(1.0,anrmx);
                    sy = amrex::Math::copysign(1.0,anrmy);
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
void mlebabeclap_ebflux (int i, int j, int k, int n,
                         Array4<Real> const& feb,
                         Array4<Real const> const& x,
                         Array4<EBCellFlag const> const& flag,
                         Array4<Real const> const& vfrc,
                         Array4<Real const> const& apx,
                         Array4<Real const> const& apy,
                         Array4<Real const> const& bc,
                         Array4<Real const> const& beb,
                         Array4<Real const> const& phieb,
                         bool is_inhomog,
                         GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real dhx = dxinv[0];

    if (!flag(i,j,k).isSingleValued())
    {
        feb(i,j,k,n) = 0.0;
    }
    else
    {
        Real kappa = vfrc(i,j,k);
        Real apxm = apx(i,j,k);
        Real apxp = apx(i+1,j,k);
        Real apym = apy(i,j,k);
        Real apyp = apy(i,j+1,k);

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
        if (amrex::Math::abs(anrmx) > amrex::Math::abs(anrmy)) {
            dg = dx_eb / amrex::Math::abs(anrmx);
            gx = bctx - dg*anrmx;
            gy = bcty - dg*anrmy;
            sx = amrex::Math::copysign(1.0,anrmx);
            sy = amrex::Math::copysign(1.0,anrmy);
        } else {
            dg = dx_eb / amrex::Math::abs(anrmy);
            gx = bctx - dg*anrmx;
            gy = bcty - dg*anrmy;
            sx = amrex::Math::copysign(1.0,anrmx);
            sy = amrex::Math::copysign(1.0,anrmy);
        }

        int ii = i - static_cast<int>(sx);
        int jj = j - static_cast<int>(sy);

        Real phig = (1.0 + gx*sx + gy*sy + gx*gy*sx*sy) * x(i ,j ,k,n)
            +       (    - gx*sx         - gx*gy*sx*sy) * x(ii,j ,k,n)
            +       (            - gy*sy - gx*gy*sx*sy) * x(i ,jj,k,n)
            +       (                    + gx*gy*sx*sy) * x(ii,jj,k,n) ;

        Real dphidn = dhx*(phib-phig)/dg;
        feb(i,j,k,n) = -beb(i,j,k,n) * dphidn;
    }
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
                       bool is_dirichlet, bool beta_on_centroid, bool phi_on_centroid,
                       Box const& vbox, int redblack, int ncomp) noexcept
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
                        int jj = j + static_cast<int>(amrex::Math::copysign(1.0,fcx(i,j,k)));
                        Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k))
                            ? amrex::Math::abs(fcx(i,j,k)) : 0.0;
                        if (!beta_on_centroid and !phi_on_centroid)
                        {
                            fxm = (1.0-fracy)*fxm + 
                                       fracy *bX(i,jj,k,n)*(phi(i,jj,k,n)-phi(i-1,jj,k,n));
                        } 
                        else if (beta_on_centroid and !phi_on_centroid)
                        {
                            fxm = (1.0-fracy)*(             -phi(i-1,j,k,n)) +
                                       fracy *(phi(i,jj,k,n)-phi(i-1,jj,k,n));
                            fxm *= bX(i,j,k,n);
                        }
                        oxm = 0.0;
                        sxm = (1.0-fracy)*sxm;
                    }

                    Real fxp =  bX(i+1,j,k,n)*phi(i+1,j,k,n);
                    Real oxp =  bX(i+1,j,k,n)*cf2;
                    Real sxp = -bX(i+1,j,k,n);
                    if (apxp != 0.0 and apxp != 1.0) {
                        int jj = j + static_cast<int>(amrex::Math::copysign(1.0,fcx(i+1,j,k)));
                        Real fracy = (ccm(i,jj,k) || ccm(i+1,jj,k))
                            ? amrex::Math::abs(fcx(i+1,j,k)) : 0.0;
                        if (!beta_on_centroid and !phi_on_centroid)
                        {
                            fxp = (1.0-fracy)*fxp + 
                                       fracy *bX(i+1,jj,k,n)*(phi(i+1,jj,k,n)-phi(i,jj,k,n));
                        } 
                        else if (beta_on_centroid and !phi_on_centroid)
                        {
                            fxp = (1.0-fracy)*(phi(i+1,j,k,n)               ) +
                                       fracy *(phi(i+1,jj,k,n)-phi(i,jj,k,n));
                            fxp *= bX(i+1,j,k,n);
                        }
                        oxp = 0.0;
                        sxp = (1.0-fracy)*sxp;
                    }

                    Real fym = -bY(i,j,k,n)*phi(i,j-1,k,n);
                    Real oym = -bY(i,j,k,n)*cf1;
                    Real sym =  bY(i,j,k,n);
                    if (apym != 0.0 and apym != 1.0) {
                        int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j,k)));
                        Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k))
                            ? amrex::Math::abs(fcy(i,j,k)) : 0.0;
                        if (!beta_on_centroid and !phi_on_centroid)
                        {
                            fym = (1.0-fracx)*fym + 
                                       fracx *bY(ii,j,k,n)*(phi(ii,j,k,n)-phi(ii,j-1,k,n));
                        } 
                        else if (beta_on_centroid and !phi_on_centroid)
                        {
                            fym = (1.0-fracx)*(             -phi( i,j-1,k,n)) +
                                       fracx *(phi(ii,j,k,n)-phi(ii,j-1,k,n));
                            fym *= bY(i,j,k,n);
                        }

                        oym = 0.0;
                        sym = (1.0-fracx)*sym;
                    }

                    Real fyp =  bY(i,j+1,k,n)*phi(i,j+1,k,n);
                    Real oyp =  bY(i,j+1,k,n)*cf3;
                    Real syp = -bY(i,j+1,k,n);
                    if (apyp != 0.0 and apyp != 1.0) {
                        int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j+1,k)));
                        Real fracx = (ccm(ii,j,k) || ccm(ii,j+1,k))
                            ? amrex::Math::abs(fcy(i,j+1,k)) : 0.0;
                        if (!beta_on_centroid and !phi_on_centroid)
                        {
                            fyp = (1.0-fracx)*fyp + 
                                       fracx*bY(ii,j+1,k,n)*(phi(ii,j+1,k,n)-phi(ii,j,k,n));
                        } 
                        else if (beta_on_centroid and !phi_on_centroid)
                        {
                            fyp = (1.0-fracx)*(phi(i,j+1,k,n)               )+
                                       fracx *(phi(ii,j+1,k,n)-phi(ii,j,k,n));
                            fyp *= bY(i,j+1,k,n);
                        }
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
                        if (amrex::Math::abs(anrmx) > amrex::Math::abs(anrmy)) {
                            dg = dx_eb / amrex::Math::abs(anrmx);
                            gx = bctx - dg*anrmx;
                            gy = bcty - dg*anrmy;
                            sx = amrex::Math::copysign(1.0,anrmx);
                            sy = amrex::Math::copysign(1.0,anrmy);
                        } else {
                            dg = dx_eb / amrex::Math::abs(anrmy);
                            gx = bctx - dg*anrmx;
                            gy = bcty - dg*anrmy;
                            sx = amrex::Math::copysign(1.0,anrmx);
                            sy = amrex::Math::copysign(1.0,anrmy);
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_flux_x (Box const& box, Array4<Real> const& fx, Array4<Real const> const& apx,
                         Array4<Real const> const& fcx, Array4<Real const> const& sol,
                         Array4<Real const> const& bX, Array4<int const> const& ccm,
                         Real dhx, int face_only, int ncomp, Box const& xbox,
                         bool beta_on_centroid, bool phi_on_centroid) noexcept
{
    int lof = xbox.smallEnd(0);
    int hif = xbox.bigEnd(0);
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (!face_only or lof == i or hif == i) {
            if (apx(i,j,k) == 0.0) {
                fx(i,j,k,n) = 0.0;
            } else if (apx(i,j,k) == 1.0) {
                fx(i,j,k,n) = -dhx*bX(i,j,k,n)*(sol(i,j,k,n)-sol(i-1,j,k,n));
            } else {
                Real fxm = bX(i,j,k,n)*(sol(i,j,k,n)-sol(i-1,j,k,n));
                int jj = j + static_cast<int>(amrex::Math::copysign(1.0,fcx(i,j,k)));
                Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? amrex::Math::abs(fcx(i,j,k)) : 0.0;
                if (!beta_on_centroid and !phi_on_centroid)
                {
                    fxm = (1.0-fracy)*fxm + fracy*bX(i,jj,k,n)*(sol(i,jj,k,n)-sol(i-1,jj,k,n));
                } else if (beta_on_centroid and !phi_on_centroid)
                    fxm = bX(i,j,k,n) * ( (1.0-fracy)*(sol(i, j,k,n)-sol(i-1, j,k,n)) +
                                               fracy *(sol(i,jj,k,n)-sol(i-1,jj,k,n)) );
                fx(i,j,k,n) = -fxm*dhx;
            }
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_flux_y (Box const& box, Array4<Real> const& fy, Array4<Real const> const& apy,
                         Array4<Real const> const& fcy, Array4<Real const> const& sol,
                         Array4<Real const> const& bY, Array4<int const> const& ccm,
                         Real dhy, int face_only, int ncomp, Box const& ybox,
                         bool beta_on_centroid, bool phi_on_centroid) noexcept
{
    int lof = ybox.smallEnd(1);
    int hif = ybox.bigEnd(1);
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (!face_only or lof == j or hif == j) {
            if (apy(i,j,k) == 0.0) {
                fy(i,j,k,n) = 0.0;
            } else if (apy(i,j,k) == 1.0) {
                fy(i,j,k,n) = -dhy*bY(i,j,k,n)*(sol(i,j,k,n)-sol(i,j-1,k,n));
            } else {
                Real fym = bY(i,j,k,n)*(sol(i,j,k,n)-sol(i,j-1,k,n));
                int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j,k)));
                Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? amrex::Math::abs(fcy(i,j,k)) : 0.0;
                if (!beta_on_centroid and !phi_on_centroid)
                {
                    fym = (1.0-fracx)*fym + fracx*bY(ii,j,k,n)*(sol(ii,j,k,n)-sol(ii,j-1,k,n));
                } else if (beta_on_centroid and !phi_on_centroid)
                    fym = bY(i,j,k,n) * ( (1.0-fracx)*(sol( i,j,k,n)-sol( i,j-1,k,n)) +
                                               fracx *(sol(ii,j,k,n)-sol(ii,j-1,k,n)) );
                fy(i,j,k,n) = -fym*dhy;
            }
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_flux_x_0 (Box const& box, Array4<Real> const& fx, Array4<Real const> const& apx,
                           Array4<Real const> const& sol, Array4<Real const> const& bX,
                           Real dhx, int face_only, int ncomp, Box const& xbox) noexcept
{
    int lof = xbox.smallEnd(0);
    int hif = xbox.bigEnd(0);
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (!face_only or lof == i or hif == i) {
            if (apx(i,j,k) == 0.0) {
                fx(i,j,k,n) = 0.0;
            } else {
                fx(i,j,k,n) = -dhx*bX(i,j,k,n)*(sol(i,j,k,n)-sol(i-1,j,k,n));
            }
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_flux_y_0 (Box const& box, Array4<Real> const& fy, Array4<Real const> const& apy,
                           Array4<Real const> const& sol, Array4<Real const> const& bY,
                           Real dhy, int face_only, int ncomp, Box const& ybox) noexcept
{
    int lof = ybox.smallEnd(1);
    int hif = ybox.bigEnd(1);
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (!face_only or lof == j or hif == j) {
            if (apy(i,j,k) == 0.0) {
                fy(i,j,k,n) = 0.0;
            } else {
                fy(i,j,k,n) = -dhy*bY(i,j,k,n)*(sol(i,j,k,n)-sol(i,j-1,k,n));
            }
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_grad_x (Box const& box, Array4<Real> const& gx, Array4<Real const> const& sol,
                         Array4<Real const> const& apx, Array4<Real const> const& fcx,
                         Array4<int const> const& ccm, 
                         Real dxi, int ncomp, bool phi_on_centroid) noexcept
{
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apx(i,j,k) == 0.0) {
            gx(i,j,k,n) = 0.0;
        } else if (apx(i,j,k) == 1.0) {
            gx(i,j,k,n) = dxi*(sol(i,j,k,n)-sol(i-1,j,k,n));
        } else {
            Real gxm = (sol(i,j,k,n)-sol(i-1,j,k,n));
            int jj = j + static_cast<int>(amrex::Math::copysign(1.0,fcx(i,j,k)));
            Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? amrex::Math::abs(fcx(i,j,k)) : 0.0;
            if (!phi_on_centroid) 
               gxm = (1.0-fracy)*gxm + fracy*(sol(i,jj,k,n)-sol(i-1,jj,k,n));
            gx(i,j,k,n) = gxm*dxi;
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_grad_y (Box const& box, Array4<Real> const& gy, Array4<Real const> const& sol,
                         Array4<Real const> const& apy, Array4<Real const> const& fcy,
                         Array4<int const> const& ccm,
                         Real dyi, int ncomp, bool phi_on_centroid) noexcept
{
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apy(i,j,k) == 0.0) {
            gy(i,j,k,n) = 0.0;
        } else if (apy(i,j,k) == 1.0) {
            gy(i,j,k,n) = dyi*(sol(i,j,k,n)-sol(i,j-1,k,n));
        } else {
            Real gym = (sol(i,j,k,n)-sol(i,j-1,k,n));
            int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j,k)));
            Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? amrex::Math::abs(fcy(i,j,k)) : 0.0;
            if (!phi_on_centroid) 
                gym = (1.0-fracx)*gym + fracx*(sol(ii,j,k,n)-sol(ii,j-1,k,n));
            gy(i,j,k,n) = gym*dyi;
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_grad_x_0 (Box const& box, Array4<Real> const& gx, Array4<Real const> const& sol,
                           Array4<Real const> const& apx, Real dxi, int ncomp) noexcept
{
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apx(i,j,k) == 0.0) {
            gx(i,j,k,n) = 0.0;
        } else {
            gx(i,j,k,n) = dxi*(sol(i,j,k,n)-sol(i-1,j,k,n));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_grad_y_0 (Box const& box, Array4<Real> const& gy, Array4<Real const> const& sol,
                           Array4<Real const> const& apy, Real dyi, int ncomp) noexcept
{
    amrex::LoopConcurrent(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apy(i,j,k) == 0.0) {
            gy(i,j,k,n) = 0.0;
        } else {
            gy(i,j,k,n) = dyi*(sol(i,j,k,n)-sol(i,j-1,k,n));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_normalize (Box const& box, Array4<Real> const& phi,
                            Real alpha, Array4<Real const> const& a,
                            Real dhx, Real dhy,
                            Array4<Real const> const& bX, Array4<Real const> const& bY,
                            Array4<const int> const& ccm, Array4<EBCellFlag const> const& flag,
                            Array4<Real const> const& vfrc,
                            Array4<Real const> const& apx, Array4<Real const> const& apy,
                            Array4<Real const> const& fcx, Array4<Real const> const& fcy,
                            Array4<Real const> const& ba, Array4<Real const> const& bc,
                            Array4<Real const> const& beb,
                            bool is_dirichlet, bool beta_on_centroid, int ncomp) noexcept
{
    amrex::Loop(box, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isRegular())
        {
            phi(i,j,k,n) /= alpha*a(i,j,k) + dhx*(bX(i,j,k,n) + bX(i+1,j,k,n))
                                           + dhy*(bY(i,j,k,n) + bY(i,j+1,k,n));
        }
        else if (flag(i,j,k).isSingleValued())
        {
            Real kappa = vfrc(i,j,k);
            Real apxm = apx(i,j,k);
            Real apxp = apx(i+1,j,k);
            Real apym = apy(i,j,k);
            Real apyp = apy(i,j+1,k);

            Real sxm =  bX(i,j,k,n);
            if (apxm != 0.0 and apxm != 1.0 && !beta_on_centroid) {
                int jj = j + static_cast<int>(amrex::Math::copysign(1.0,fcx(i,j,k)));
                Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k))
                    ? amrex::Math::abs(fcx(i,j,k)) : 0.0;
                sxm = (1.0-fracy)*sxm;
            }

            Real sxp = -bX(i+1,j,k,n);
            if (apxp != 0.0 and apxp != 1.0 && !beta_on_centroid) {
                int jj = j + static_cast<int>(amrex::Math::copysign(1.0,fcx(i+1,j,k)));
                Real fracy = (ccm(i,jj,k) || ccm(i+1,jj,k))
                    ? amrex::Math::abs(fcx(i+1,j,k)) : 0.0;
                sxp = (1.0-fracy)*sxp;
            }

            Real sym =  bY(i,j,k,n);
            if (apym != 0.0 and apym != 1.0 && !beta_on_centroid) {
                int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j,k)));
                Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k))
                    ? amrex::Math::abs(fcy(i,j,k)) : 0.0;
                sym = (1.0-fracx)*sym;
            }

            Real syp = -bY(i,j+1,k,n);
            if (apyp != 0.0 and apyp != 1.0 && !beta_on_centroid) {
                int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j+1,k)));
                Real fracx = (ccm(ii,j,k) || ccm(ii,j+1,k))
                    ? amrex::Math::abs(fcy(i,j+1,k)) : 0.0;
                syp = (1.0-fracx)*syp;
            }

            Real vfrcinv = (1.0/kappa);
            Real gamma = alpha*a(i,j,k) + vfrcinv *
                (dhx*(apxm*sxm-apxp*sxp) +
                 dhy*(apym*sym-apyp*syp));

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
                if (amrex::Math::abs(anrmx) > amrex::Math::abs(anrmy)) {
                    dg = dx_eb / amrex::Math::abs(anrmx);
                    gx = bctx - dg*anrmx;
                    gy = bcty - dg*anrmy;
                    sx = amrex::Math::copysign(1.0,anrmx);
                    sy = amrex::Math::copysign(1.0,anrmy);
                } else {
                    dg = dx_eb / amrex::Math::abs(anrmy);
                    gx = bctx - dg*anrmx;
                    gy = bcty - dg*anrmy;
                    sx = amrex::Math::copysign(1.0,anrmx);
                    sy = amrex::Math::copysign(1.0,anrmy);
                }

                Real phig_gamma = (1.0 + gx*sx + gy*sy + gx*gy*sx*sy);
                Real feb_gamma = -phig_gamma/dg * ba(i,j,k) * beb(i,j,k,n);
                gamma += vfrcinv*(-dhx)*feb_gamma;
            }

            phi(i,j,k,n) /= gamma;
        }
    });
}

}
#endif
