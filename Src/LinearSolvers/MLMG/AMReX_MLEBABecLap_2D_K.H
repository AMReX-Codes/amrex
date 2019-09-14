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
                Real dx_eb = amrex::max(0.3,(kappa*kappa-0.25)/(2.0*kappa));

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

}

#endif
