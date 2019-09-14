#ifndef AMREX_MLEBABECLAP_3D_K_H_
#define AMREX_MLEBABECLAP_3D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebabeclap_adotx (Box const& box, Array4<Real> const& y,
                        Array4<Real const> const& x, Array4<Real const> const& a,
                        Array4<Real const> const& bX, Array4<Real const> const& bY,
                        Array4<Real const> const& bZ, Array4<const int> const& ccm,
                        Array4<EBCellFlag const> const& flag,
                        Array4<Real const> const& vfrc, Array4<Real const> const& apx,
                        Array4<Real const> const& apy, Array4<Real const> const& apz,
                        Array4<Real const> const& fcx, Array4<Real const> const& fcy,
                        Array4<Real const> const& fcz, Array4<Real const> const& ba,
                        Array4<Real const> const& bc, Array4<Real const> const& beb,
                        bool is_dirichlet, Array4<Real const> const& phieb,
                        bool is_inhomog, GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                        Real alpha, Real beta, int ncomp)
{
    Real dhx = beta*dxinv[0]*dxinv[0];
    Real dhy = beta*dxinv[1]*dxinv[1];
    Real dhz = beta*dxinv[2]*dxinv[2];

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
                        -bX(i  ,j,k,n)*(x(i  ,j,k,n) - x(i-1,j,k,n)))
                - dhy * (bY(i,j+1,k,n)*(x(i,j+1,k,n) - x(i,j  ,k,n))
                        -bY(i,j  ,k,n)*(x(i,j  ,k,n) - x(i,j-1,k,n)))
                - dhz * (bZ(i,j,k+1,n)*(x(i,j,k+1,n) - x(i,j,k  ,n))
                        -bZ(i,j,k  ,n)*(x(i,j,k  ,n) - x(i,j,k-1,n)));
        }
        else
        {
            Real kappa = vfrc(i,j,k);
            Real apxm = apx(i,j,k);
            Real apxp = apx(i+1,j,k);
            Real apym = apy(i,j,k);
            Real apyp = apy(i,j+1,k);
            Real apzm = apz(i,j,k);
            Real apzp = apz(i,j,k+1);

            Real fxm = bX(i,j,k,n)*(x(i,j,k,n) - x(i-1,j,k,n));
            if (apxm != 0.0 and apxm != 1.0) {
                int jj = j + static_cast<int>(std::copysign(1.0, fcx(i,j,k,0)));
                int kk = k + static_cast<int>(std::copysign(1.0, fcx(i,j,k,1)));
                Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? std::abs(fcx(i,j,k,0)) : 0.0;
                Real fracz = (ccm(i-1,j,kk) || ccm(i,j,kk)) ? std::abs(fcx(i,j,k,1)) : 0.0;
                fxm = (1.0-fracy)*(1.0-fracz)*fxm +
                    fracy*(1.0-fracz)*bX(i,jj,k ,n)*(x(i,jj,k ,n)-x(i-1,jj,k ,n)) +
                    fracz*(1.0-fracy)*bX(i,j ,kk,n)*(x(i,j ,kk,n)-x(i-1,j ,kk,n)) +
                    fracy*     fracz *bX(i,jj,kk,n)*(x(i,jj,kk,n)-x(i-1,jj,kk,n));
            }

            Real fxp = bX(i+1,j,k,n)*(x(i+1,j,k,n) - x(i,j,k,n));
            if (apxp != 0.0 and apxp != 1.0) {
                int jj = j + static_cast<int>(std::copysign(1.0,fcx(i+1,j,k,0)));
                int kk = k + static_cast<int>(std::copysign(1.0,fcx(i+1,j,k,1)));
                Real fracy = (ccm(i,jj,k) || ccm(i+1,jj,k)) ? std::abs(fcx(i+1,j,k,0)) : 0.0;
                Real fracz = (ccm(i,j,kk) || ccm(i+1,j,kk)) ? std::abs(fcx(i+1,j,k,1)) : 0.0;
                fxp = (1.0-fracy)*(1.0-fracz)*fxp +
                    fracy*(1.0-fracz)*bX(i+1,jj,k ,n)*(x(i+1,jj,k ,n)-x(i,jj,k ,n)) +
                    fracz*(1.0-fracy)*bX(i+1,j ,kk,n)*(x(i+1,j ,kk,n)-x(i,j ,kk,n)) +
                    fracy*     fracz *bX(i+1,jj,kk,n)*(x(i+1,jj,kk,n)-x(i,jj,kk,n));
            }

            Real fym = bY(i,j,k,n)*(x(i,j,k,n) - x(i,j-1,k,n));
            if (apym != 0.0 and apym != 1.0) {
                int ii = i + static_cast<int>(std::copysign(1.0,fcy(i,j,k,0)));
                int kk = k + static_cast<int>(std::copysign(1.0,fcy(i,j,k,1)));
                Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? std::abs(fcy(i,j,k,0)) : 0.0;
                Real fracz = (ccm(i,j-1,kk) || ccm(i,j,kk)) ? std::abs(fcy(i,j,k,1)) : 0.0;
                fym = (1.0-fracx)*(1.0-fracz)*fym +
                    fracx*(1.0-fracz)*bY(ii,j,k ,n)*(x(ii,j,k ,n)-x(ii,j-1,k ,n)) +
                    fracz*(1.0-fracx)*bY(i ,j,kk,n)*(x(i ,j,kk,n)-x(i ,j-1,kk,n)) +
                    fracx*     fracz *bY(ii,j,kk,n)*(x(ii,j,kk,n)-x(ii,j-1,kk,n));
            }

            Real fyp = bY(i,j+1,k,n)*(x(i,j+1,k,n) - x(i,j,k,n));
            if (apyp != 0.0 and apyp != 1.0) {
                int ii = i + static_cast<int>(std::copysign(1.0,fcy(i,j+1,k,0)));
                int kk = k + static_cast<int>(std::copysign(1.0,fcy(i,j+1,k,1)));
                Real fracx = (ccm(ii,j,k) || ccm(ii,j+1,k)) ? std::abs(fcy(i,j+1,k,0)) : 0.0;
                Real fracz = (ccm(i,j,kk) || ccm(i,j+1,kk)) ? std::abs(fcy(i,j+1,k,1)) : 0.0;
                fyp = (1.0-fracx)*(1.0-fracz)*fyp +
                    fracx*(1.0-fracz)*bY(ii,j+1,k ,n)*(x(ii,j+1,k ,n)-x(ii,j,k ,n)) +
                    fracz*(1.0-fracx)*bY(i ,j+1,kk,n)*(x(i ,j+1,kk,n)-x(i ,j,kk,n)) +
                    fracx*     fracz *bY(ii,j+1,kk,n)*(x(ii,j+1,kk,n)-x(ii,j,kk,n));
            }

            Real fzm = bZ(i,j,k,n)*(x(i,j,k,n) - x(i,j,k-1,n));
            if (apzm != 0.0 and apzm != 1.0) {
                int ii = i + static_cast<int>(std::copysign(1.0,fcz(i,j,k,0)));
                int jj = j + static_cast<int>(std::copysign(1.0,fcz(i,j,k,1)));
                Real fracx = (ccm(ii,j,k-1) || ccm(ii,j,k)) ? std::abs(fcz(i,j,k,0)) : 0.0;
                Real fracy = (ccm(i,jj,k-1) || ccm(i,jj,k)) ? std::abs(fcz(i,j,k,1)) : 0.0;
                fzm = (1.0-fracx)*(1.0-fracy)*fzm +
                    fracx*(1.0-fracy)*bZ(ii,j ,k,n)*(x(ii,j ,k,n)-x(ii,j ,k-1,n)) +
                    fracy*(1.0-fracx)*bZ(i ,jj,k,n)*(x(i ,jj,k,n)-x(i ,jj,k-1,n)) +
                    fracx*     fracy *bZ(ii,jj,k,n)*(x(ii,jj,k,n)-x(ii,jj,k-1,n));
            }

            Real fzp = bZ(i,j,k+1,n)*(x(i,j,k+1,n) - x(i,j,k,n));
            if (apzp != 0.0 and apzp != 1.0) {
                int ii = i + static_cast<int>(std::copysign(1.0,fcz(i,j,k+1,0)));
                int jj = j + static_cast<int>(std::copysign(1.0,fcz(i,j,k+1,1)));
                Real fracx = (ccm(ii,j,k) || ccm(ii,j,k+1)) ? std::abs(fcz(i,j,k+1,0)) : 0.0;
                Real fracy = (ccm(i,jj,k) || ccm(i,jj,k+1)) ? std::abs(fcz(i,j,k+1,1)) : 0.0;
                fzp = (1.0-fracx)*(1.0-fracy)*fzp +
                    fracx*(1.0-fracy)*bZ(ii,j ,k+1,n)*(x(ii,j ,k+1,n)-x(ii,j ,k,n)) +
                    fracy*(1.0-fracx)*bZ(i ,jj,k+1,n)*(x(i ,jj,k+1,n)-x(i ,jj,k,n)) +
                    fracx*     fracy *bZ(ii,jj,k+1,n)*(x(ii,jj,k+1,n)-x(ii,jj,k,n));
            }

            Real feb = 0.0;
            if (is_dirichlet) {
                Real dapx = apxm-apxp;
                Real dapy = apym-apyp;
                Real dapz = apzm-apzp;
                Real anorm = std::sqrt(dapx*dapx+dapy*dapy+dapz+dapz);
                Real anorminv = 1.0/anorm;
                Real anrmx = dapx * anorminv;
                Real anrmy = dapy * anorminv;
                Real anrmz = dapz * anorminv;

                Real phib = is_inhomog ? phieb(i,j,k,n) : 0.0;

                Real bctx = bc(i,j,k,0);
                Real bcty = bc(i,j,k,1);
                Real bctz = bc(i,j,k,2);
                Real dx_eb = amrex::max(0.3,(kappa*kappa-0.25)/(2.0*kappa));

                Real dg = dx_eb / amrex::max(std::abs(anrmx), std::abs(anrmy),
                                             std::abs(anrmz));
                Real gx = bctx - dg*anrmx;
                Real gy = bcty - dg*anrmy;
                Real gz = bctz - dg*anrmz;
                Real sx = std::copysign(1.0,anrmx);
                Real sy = std::copysign(1.0,anrmy);
                Real sz = std::copysign(1.0,anrmz);
                int ii = i - static_cast<int>(sx);
                int jj = j - static_cast<int>(sy);
                int kk = k - static_cast<int>(sz);

                gx = sx*gx;
                gy = sy*gy;
                gz = sz*gz;
                Real gxy = gx*gy;
                Real gxz = gx*gz;
                Real gyz = gy*gz;
                Real gxyz = gx*gy*gz;
                Real phig = (1.0+gx+gy+gz+gxy+gxz+gyz+gxyz) * x(i ,j ,k ,n)
                    +       (-gz - gxz - gyz - gxyz)        * x(i ,j ,kk,n)
                    +       (-gy - gxy - gyz - gxyz)        * x(i ,jj,k ,n)
                    +       (gyz + gxyz)                    * x(i ,jj,kk,n)
                    +       (-gx - gxy - gxz - gxyz)        * x(ii,j ,k ,n)
                    +       (gxz + gxyz)                    * x(ii,j ,kk,n)
                    +       (gxy + gxyz)                    * x(ii,jj,k ,n)
                    +       (-gxyz)                         * x(ii,jj,kk,n);

                Real dphidn = (phib-phig)/dg;

                feb = dphidn * ba(i,j,k) * beb(i,j,k,n);
            }

            y(i,j,k,n) = alpha*a(i,j,k)*x(i,j,k,n) + (1.0/kappa) *
                (dhx*(apxm*fxm - apxp*fxp) +
                 dhy*(apym*fym - apyp*fyp) +
                 dhz*(apzm*fzm - apzp*fzp) - dhx*feb);
        }
    });
}

}

#endif
