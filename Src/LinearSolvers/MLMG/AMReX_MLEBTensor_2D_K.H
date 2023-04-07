#ifndef AMREX_ML_EB_TENSOR_2D_K_H_
#define AMREX_ML_EB_TENSOR_2D_K_H_
#include <AMReX_Config.H>

#include <AMReX_MLEBABecLap_K.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms_fx (Box const& box, Array4<Real> const& fx,
                                Array4<Real const> const& vel,
                                Array4<Real const> const& etax,
                                Array4<Real const> const& kapx,
                                Array4<Real const> const& apx,
                                Array4<EBCellFlag const> const& flag,
                                GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apx(i,j,0) == Real(0.0))
            {
                fx(i,j,0,0) = Real(0.0);
                fx(i,j,0,1) = Real(0.0);
            }
            else
            {
                int jhip = j + flag(i  ,j,0).isConnected(0, 1,0);
                int jhim = j - flag(i  ,j,0).isConnected(0,-1,0);
                int jlop = j + flag(i-1,j,0).isConnected(0, 1,0);
                int jlom = j - flag(i-1,j,0).isConnected(0,-1,0);
                Real whi = mlebtensor_weight(jhip-jhim);
                Real wlo = mlebtensor_weight(jlop-jlom);
                Real dudy = mlebtensor_dy_on_xface(i,j,k,0,vel,dyi,
                                                   whi,wlo,jhip,jhim,jlop,jlom);
                Real dvdy = mlebtensor_dy_on_xface(i,j,k,1,vel,dyi,
                                                   whi,wlo,jhip,jhim,jlop,jlom);
                Real divu = dvdy;
                Real xif = kapx(i,j,0);
                Real mun = Real(0.75)*(etax(i,j,0,0)-xif);// restore the original eta
                Real mut =             etax(i,j,0,1);
                fx(i,j,0,0) = -mun*(-twoThirds*divu) - xif*divu;
                fx(i,j,0,1) = -mut*dudy;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms_fy (Box const& box, Array4<Real> const& fy,
                                Array4<Real const> const& vel,
                                Array4<Real const> const& etay,
                                Array4<Real const> const& kapy,
                                Array4<Real const> const& apy,
                                Array4<EBCellFlag const> const& flag,
                                GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apy(i,j,0) == Real(0.0))
            {
                fy(i,j,0,0) = Real(0.0);
                fy(i,j,0,1) = Real(0.0);
            }
            else
            {
                int ihip = i + flag(i,j  ,0).isConnected( 1,0,0);
                int ihim = i - flag(i,j  ,0).isConnected(-1,0,0);
                int ilop = i + flag(i,j-1,0).isConnected( 1,0,0);
                int ilom = i - flag(i,j-1,0).isConnected(-1,0,0);
                Real whi = mlebtensor_weight(ihip-ihim);
                Real wlo = mlebtensor_weight(ilop-ilom);
                Real dudx = mlebtensor_dx_on_yface(i,j,k,0,vel,dxi,
                                                   whi,wlo,ihip,ihim,ilop,ilom);
                Real dvdx = mlebtensor_dx_on_yface(i,j,k,1,vel,dxi,
                                                   whi,wlo,ihip,ihim,ilop,ilom);
                Real divu = dudx;
                Real xif = kapy(i,j,0);
                Real mun = Real(0.75)*(etay(i,j,0,1)-xif);// restore the original eta
                Real mut =             etay(i,j,0,0);
                fy(i,j,0,0) = -mut*dvdx;
                fy(i,j,0,1) = -mun*(-twoThirds*divu) - xif*divu;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms_fx (Box const& box, Array4<Real> const& fx,
                                Array4<Real const> const& vel,
                                Array4<Real const> const& etax,
                                Array4<Real const> const& kapx,
                                Array4<Real const> const& apx,
                                Array4<EBCellFlag const> const& flag,
                                GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                Array4<Real const> const& bvxlo,
                                Array4<Real const> const& bvxhi,
                                Array2D<BoundCond,
                                        0,2*AMREX_SPACEDIM,
                                        0,AMREX_SPACEDIM> const& bct,
                                Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apx(i,j,0) == Real(0.0))
            {
                fx(i,j,0,0) = Real(0.0);
                fx(i,j,0,1) = Real(0.0);
            }
            else
            {
                int jhip = j + flag(i  ,j,0).isConnected(0, 1,0);
                int jhim = j - flag(i  ,j,0).isConnected(0,-1,0);
                int jlop = j + flag(i-1,j,0).isConnected(0, 1,0);
                int jlom = j - flag(i-1,j,0).isConnected(0,-1,0);
                Real whi = mlebtensor_weight(jhip-jhim);
                Real wlo = mlebtensor_weight(jlop-jlom);
                Real dudy = mlebtensor_dy_on_xface(i,j,k,0,vel,dyi,
                                                   bvxlo,bvxhi,bct,dlo,dhi,
                                                   whi,wlo,jhip,jhim,jlop,jlom);
                Real dvdy = mlebtensor_dy_on_xface(i,j,k,1,vel,dyi,
                                                   bvxlo,bvxhi,bct,dlo,dhi,
                                                   whi,wlo,jhip,jhim,jlop,jlom);
                Real divu = dvdy;
                Real xif = kapx(i,j,0);
                Real mun = Real(0.75)*(etax(i,j,0,0)-xif);// restore the original eta
                Real mut =       etax(i,j,0,1);
                fx(i,j,0,0) = -mun*(-twoThirds*divu) - xif*divu;
                fx(i,j,0,1) = -mut*dudy;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms_fy (Box const& box, Array4<Real> const& fy,
                                Array4<Real const> const& vel,
                                Array4<Real const> const& etay,
                                Array4<Real const> const& kapy,
                                Array4<Real const> const& apy,
                                Array4<EBCellFlag const> const& flag,
                                GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                Array4<Real const> const& bvylo,
                                Array4<Real const> const& bvyhi,
                                Array2D<BoundCond,
                                        0,2*AMREX_SPACEDIM,
                                        0,AMREX_SPACEDIM> const& bct,
                                Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apy(i,j,0) == Real(0.0))
            {
                fy(i,j,0,0) = Real(0.0);
                fy(i,j,0,1) = Real(0.0);
            }
            else
            {
                int ihip = i + flag(i,j  ,0).isConnected( 1,0,0);
                int ihim = i - flag(i,j  ,0).isConnected(-1,0,0);
                int ilop = i + flag(i,j-1,0).isConnected( 1,0,0);
                int ilom = i - flag(i,j-1,0).isConnected(-1,0,0);
                Real whi = mlebtensor_weight(ihip-ihim);
                Real wlo = mlebtensor_weight(ilop-ilom);
                Real dudx = mlebtensor_dx_on_yface(i,j,k,0,vel,dxi,
                                                   bvylo,bvyhi,bct,dlo,dhi,
                                                   whi,wlo,ihip,ihim,ilop,ilom);
                Real dvdx = mlebtensor_dx_on_yface(i,j,k,1,vel,dxi,
                                                   bvylo,bvyhi,bct,dlo,dhi,
                                                   whi,wlo,ihip,ihim,ilop,ilom);
                Real divu = dudx;
                Real xif = kapy(i,j,0);
                Real mun = Real(0.75)*(etay(i,j,0,1)-xif);// restore the original eta
                Real mut =             etay(i,j,0,0);
                fy(i,j,0,0) = -mut*dvdx;
                fy(i,j,0,1) = -mun*(-twoThirds*divu) - xif*divu;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms (Box const& box, Array4<Real> const& Ax,
                             Array4<Real const> const& fx,
                             Array4<Real const> const& fy,
                             Array4<Real const> const& vel,
                             Array4<Real const> const& velb,
                             Array4<Real const> const& etab,
                             Array4<Real const> const& kapb,
                             Array4<int const> const& ccm,
                             Array4<EBCellFlag const> const& flag,
                             Array4<Real const> const& vol,
                             Array4<Real const> const& apx,
                             Array4<Real const> const& apy,
                             Array4<Real const> const& fcx,
                             Array4<Real const> const& fcy,
                             Array4<Real const> const& bc,
                             bool is_dirichlet, bool is_inhomog,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                             Real bscalar) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if (flag(i,j,0).isRegular())
            {
                Ax(i,j,0,0) += bscalar*(dxi*(fx(i+1,j  ,0,0) - fx(i,j,0,0))
                                      + dyi*(fy(i  ,j+1,0,0) - fy(i,j,0,0)));
                Ax(i,j,0,1) += bscalar*(dxi*(fx(i+1,j  ,0,1) - fx(i,j,0,1))
                                      + dyi*(fy(i  ,j+1,0,1) - fy(i,j,0,1)));
            }
            else if (flag(i,j,0).isSingleValued())
            {
                Real fxm_0 = fx(i,j,0,0);
                Real fxm_1 = fx(i,j,0,1);
                if (apx(i,j,0) > Real(0.0) && apx(i,j,0) < Real(1.0)) {
                    int jj = j + (fcx(i,j,0,0) >= Real(0.0) ? 1 : -1);
                    Real fracy = (ccm(i-1,jj,0) || ccm(i,jj,0)) ? std::abs(fcx(i,j,0,0)) : Real(0.0);
                    fxm_0 = (Real(1.0)-fracy)*fxm_0 + fracy*fx(i,jj,0,0);
                    fxm_1 = (Real(1.0)-fracy)*fxm_1 + fracy*fx(i,jj,0,1);
                }

                Real fxp_0 = fx(i+1,j,0,0);
                Real fxp_1 = fx(i+1,j,0,1);
                if (apx(i+1,j,0) > Real(0.0) && apx(i+1,j,0) < Real(1.0)) {
                    int jj = j + (fcx(i+1,j,0,0) >= Real(0.0) ? 1 : -1);
                    Real fracy = (ccm(i,jj,0) || ccm(i+1,jj,0)) ? std::abs(fcx(i+1,j,0,0)) : Real(0.0);
                    fxp_0 = (Real(1.0)-fracy)*fxp_0 + fracy*fx(i+1,jj,0,0);
                    fxp_1 = (Real(1.0)-fracy)*fxp_1 + fracy*fx(i+1,jj,0,1);
                }

                Real fym_0 = fy(i,j,0,0);
                Real fym_1 = fy(i,j,0,1);
                if (apy(i,j,0) > Real(0.0) && apy(i,j,0) < Real(1.0)) {
                    int ii = i + (fcy(i,j,0,0) >= Real(0.0) ? 1 : -1);
                    Real fracx = (ccm(ii,j-1,0) || ccm(ii,j,0)) ? std::abs(fcy(i,j,0,0)) : Real(0.0);
                    fym_0 = (Real(1.0)-fracx)*fym_0 + fracx*fy(ii,j,0,0);
                    fym_1 = (Real(1.0)-fracx)*fym_1 + fracx*fy(ii,j,0,1);
                }

                Real fyp_0 = fy(i,j+1,0,0);
                Real fyp_1 = fy(i,j+1,0,1);
                if (apy(i,j+1,0) > Real(0.0) && apy(i,j+1,0) < Real(1.0)) {
                    int ii = i + (fcy(i,j+1,0,0) >= Real(0.0) ? 1 : -1);
                    Real fracx = (ccm(ii,j,0) || ccm(ii,j+1,0)) ? std::abs(fcy(i,j+1,0,0)) : Real(0.0);
                    fyp_0 = (Real(1.0)-fracx)*fyp_0 + fracx*fy(ii,j+1,0,0);
                    fyp_1 = (Real(1.0)-fracx)*fyp_1 + fracx*fy(ii,j+1,0,1);
                }

                Real kappa = vol(i,j,0);
                Real feb_0 = Real(0.0);
                Real feb_1 = Real(0.0);
                if (is_dirichlet)
                {
                    Real dapx = apx(i+1,j,0)-apx(i,j,0);
                    Real dapy = apy(i,j+1,0)-apy(i,j,0);
                    Real anorminv = Real(1.0)/std::sqrt(dapx*dapx+dapy*dapy);
                    Real anrmx = -dapx * anorminv;
                    Real anrmy = -dapy * anorminv;

                    Real velb_0 = 0, velb_1 = 0;

                    if (is_inhomog) {
                        velb_0 = velb(i,j,0,0);
                        velb_1 = velb(i,j,0,1);
                    }

                    Real dx_eb = amrex::max(Real(0.3), (kappa*kappa-Real(0.25))/(Real(2.)*kappa));
                    Real dg = dx_eb / amrex::max(std::abs(anrmx),std::abs(anrmy));
                    Real dginv = Real(1.0)/dg;
                    Real gx = bc(i,j,0,0) - dg*anrmx;
                    Real gy = bc(i,j,0,1) - dg*anrmy;
                    int isx = (anrmx > Real(0.0)) ? 1 : -1;
                    int isy = (anrmy > Real(0.0)) ? 1 : -1;
                    int ii = i - isx;
                    int jj = j - isy;

                    gx *= isx;
                    gy *= isy;
                    Real gxy = gx*gy;
                    Real oneggg = Real(1.0)+gx+gy+gxy;

                    Real velg = oneggg * vel(i ,j ,0,0)
                        +  (-gy - gxy) * vel(i ,jj,0,0)
                        +  (-gx - gxy) * vel(ii,j ,0,0)
                        +    gxy       * vel(ii,jj,0,0);
                    Real dudn = (velb_0-velg) * dginv;

                    velg      = oneggg * vel(i ,j ,0,1)
                        +  (-gy - gxy) * vel(i ,jj,0,1)
                        +  (-gx - gxy) * vel(ii,j ,0,1)
                        +    gxy       * vel(ii,jj,0,1);
                    Real dvdn = (velb_1-velg) * dginv;

                    // transverse derivatives are zero on the wall
                    Real dudx = dudn * anrmx;
                    Real dudy = dudn * anrmy;
                    Real dvdx = dvdn * anrmx;
                    Real dvdy = dvdn * anrmy;
                    Real divu = dudx + dvdy;
                    Real xi = kapb(i,j,0);
                    Real mu = etab(i,j,0);
                    Real tautmp = (xi-(2./3.)*mu)*divu;
                    // Note that mu*(grad v) has been included already in MLEBABecLap
                    Real tauxx = mu*dudx + tautmp;
                    Real tauyx = mu*dvdx;
                    Real tauyy = mu*dvdy + tautmp;
                    Real tauxy = mu*dudy;
                    // assuming dxi == dyi == dzi
                    feb_0 = dxi*(dapx*tauxx + dapy*tauyx);
                    feb_1 = dxi*(dapx*tauxy + dapy*tauyy);
                }

                Real volinv = bscalar / kappa;
                Ax(i,j,0,0) += volinv * (dxi*(apx(i,j,0)*fxm_0-apx(i+1,j,0)*fxp_0)
                                        +dyi*(apy(i,j,0)*fym_0-apy(i,j+1,0)*fyp_0)
                                        -dxi*                               feb_0);
                Ax(i,j,0,1) += volinv * (dxi*(apx(i,j,0)*fxm_1-apx(i+1,j,0)*fxp_1)
                                        +dyi*(apy(i,j,0)*fym_1-apy(i,j+1,0)*fyp_1)
                                        -dxi*                               feb_1);
            }
        }
    }
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_flux_0 (Box const& box,
                        Array4<Real> const& Ax,
                        Array4<Real const> const& fx,
                        Array4<Real const> const& ap,
                        Real bscalar) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
          if (ap(i,j,k) != Real(0.0)) {
              for (int n=0; n<AMREX_SPACEDIM; n++)
                  Ax(i,j,k,n) += bscalar*fx(i,j,k,n);
          }
          //else MLEBABec::compFlux already set Ax = 0
        }
    }
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_flux_x (Box const& box, Array4<Real> const& Ax,
                        Array4<Real const> const& fx, Array4<Real const> const& apx,
                        Array4<Real const> const& fcx,
                        Real const bscalar, Array4<int const> const& ccm,
                        int face_only, Box const& xbox) noexcept
{
    int lof = xbox.smallEnd(0);
    int hif = xbox.bigEnd(0);
    amrex::LoopConcurrent(box, [=] (int i, int j, int k) noexcept
    {
        if (!face_only || lof == i || hif == i) {
          if (apx(i,j,k) == Real(1.0)) {
            for (int n=0; n<AMREX_SPACEDIM; n++)
                Ax(i,j,k,n) += bscalar*fx(i,j,k,n);
          }
          else if (apx(i,j,k) != 0.) {
                Real fxm_0 = fx(i,j,0,0);
                Real fxm_1 = fx(i,j,0,1);

                int jj = j + (fcx(i,j,0,0) >= Real(0.0) ? 1 : -1);
                Real fracy = (ccm(i-1,jj,0) || ccm(i,jj,0)) ? std::abs(fcx(i,j,0,0)) : Real(0.0);
                fxm_0 = (Real(1.0)-fracy)*fxm_0 + fracy*fx(i,jj,0,0);
                fxm_1 = (Real(1.0)-fracy)*fxm_1 + fracy*fx(i,jj,0,1);

                Ax(i,j,k,0) += bscalar*fxm_0;
                Ax(i,j,k,1) += bscalar*fxm_1;
          }
          //else MLEBABec::compFlux already set Ax = 0
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_flux_y (Box const& box, Array4<Real> const& Ay,
                        Array4<Real const> const& fy, Array4<Real const> const& apy,
                        Array4<Real const> const& fcy,
                        Real const bscalar, Array4<int const> const& ccm,
                        int face_only, Box const& ybox) noexcept
{
    int lof = ybox.smallEnd(1);
    int hif = ybox.bigEnd(1);
    amrex::LoopConcurrent(box, [=] (int i, int j, int k) noexcept
    {
        if (!face_only || lof == j || hif == j) {
          if (apy(i,j,k) == Real(1.0)) {
            for (int n=0; n<AMREX_SPACEDIM; n++)
                Ay(i,j,k,n) += bscalar*fy(i,j,k,n);
          }
          else if (apy(i,j,k) != 0.) {
                Real fym_0 = fy(i,j,0,0);
                Real fym_1 = fy(i,j,0,1);

                int ii = i + (fcy(i,j,0,0) >= Real(0.0) ? 1 : -1);
                Real fracx = (ccm(ii,j-1,0) || ccm(ii,j,0)) ? std::abs(fcy(i,j,0,0)) : Real(0.0);
                fym_0 = (Real(1.0)-fracx)*fym_0 + fracx*fy(ii,j,0,0);
                fym_1 = (Real(1.0)-fracx)*fym_1 + fracx*fy(ii,j,0,1);

                Ay(i,j,k,0) += bscalar*fym_0;
                Ay(i,j,k,1) += bscalar*fym_1;
          }
          //else MLEBABec::compFlux already set Ay = 0
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel, Array4<Real const> const& apx,
                              Array4<EBCellFlag const> const& flag,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apx(i,j,0) == Real(0.0))
            {
                fx(i,j,0,0) = Real(0.0);
                fx(i,j,0,1) = Real(0.0);
                fx(i,j,0,2) = Real(0.0);
                fx(i,j,0,3) = Real(0.0);
            }
            else
            {
                Real dudx = (vel(i,j,0,0) - vel(i-1,j,0,0))*dxi;
                Real dvdx = (vel(i,j,0,1) - vel(i-1,j,0,1))*dxi;

                int jhip = j + flag(i  ,j,0).isConnected(0, 1,0);
                int jhim = j - flag(i  ,j,0).isConnected(0,-1,0);
                int jlop = j + flag(i-1,j,0).isConnected(0, 1,0);
                int jlom = j - flag(i-1,j,0).isConnected(0,-1,0);
                Real whi = mlebtensor_weight(jhip-jhim);
                Real wlo = mlebtensor_weight(jlop-jlom);
                Real dudy = mlebtensor_dy_on_xface(i,j,k,0,vel,dyi,
                                                   whi,wlo,jhip,jhim,jlop,jlom);
                Real dvdy = mlebtensor_dy_on_xface(i,j,k,1,vel,dyi,
                                                   whi,wlo,jhip,jhim,jlop,jlom);
                fx(i,j,0,0) = dudx;
                fx(i,j,0,1) = dvdx;
                fx(i,j,0,2) = dudy;
                fx(i,j,0,3) = dvdy;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel, Array4<Real const> const& apy,
                              Array4<EBCellFlag const> const& flag,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apy(i,j,0) == 0.) {
                fy(i,j,0,0) = Real(0.0);
                fy(i,j,0,1) = Real(0.0);
                fy(i,j,0,2) = Real(0.0);
                fy(i,j,0,3) = Real(0.0);
            }
            else
            {
                int ihip = i + flag(i,j  ,0).isConnected( 1,0,0);
                int ihim = i - flag(i,j  ,0).isConnected(-1,0,0);
                int ilop = i + flag(i,j-1,0).isConnected( 1,0,0);
                int ilom = i - flag(i,j-1,0).isConnected(-1,0,0);
                Real whi = mlebtensor_weight(ihip-ihim);
                Real wlo = mlebtensor_weight(ilop-ilom);
                Real dudx = mlebtensor_dx_on_yface(i,j,k,0,vel,dxi,
                                                   whi,wlo,ihip,ihim,ilop,ilom);
                Real dvdx = mlebtensor_dx_on_yface(i,j,k,1,vel,dxi,
                                                   whi,wlo,ihip,ihim,ilop,ilom);

                Real dudy = (vel(i,j,0,0) - vel(i,j-1,0,0))*dyi;
                Real dvdy = (vel(i,j,0,1) - vel(i,j-1,0,1))*dyi;
                fy(i,j,0,0) = dudx;
                fy(i,j,0,1) = dvdx;
                fy(i,j,0,2) = dudy;
                fy(i,j,0,3) = dvdy;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel, Array4<Real const> const& apx,
                              Array4<EBCellFlag const> const& flag,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvxlo,
                              Array4<Real const> const& bvxhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apx(i,j,0) == Real(0.0))
            {
                fx(i,j,0,0) = Real(0.0);
                fx(i,j,0,1) = Real(0.0);
                fx(i,j,0,2) = Real(0.0);
                fx(i,j,0,3) = Real(0.0);
            }
            else
            {
                Real dudx = (vel(i,j,0,0) - vel(i-1,j,0,0))*dxi;
                Real dvdx = (vel(i,j,0,1) - vel(i-1,j,0,1))*dxi;

                int jhip = j + flag(i  ,j,0).isConnected(0, 1,0);
                int jhim = j - flag(i  ,j,0).isConnected(0,-1,0);
                int jlop = j + flag(i-1,j,0).isConnected(0, 1,0);
                int jlom = j - flag(i-1,j,0).isConnected(0,-1,0);
                Real whi = mlebtensor_weight(jhip-jhim);
                Real wlo = mlebtensor_weight(jlop-jlom);
                Real dudy = mlebtensor_dy_on_xface(i,j,k,0,vel,dyi,
                                                   bvxlo,bvxhi,bct,dlo,dhi,
                                                   whi,wlo,jhip,jhim,jlop,jlom);
                Real dvdy = mlebtensor_dy_on_xface(i,j,k,1,vel,dyi,
                                                   bvxlo,bvxhi,bct,dlo,dhi,
                                                   whi,wlo,jhip,jhim,jlop,jlom);
                fx(i,j,0,0) = dudx;
                fx(i,j,0,1) = dvdx;
                fx(i,j,0,2) = dudy;
                fx(i,j,0,3) = dvdy;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel, Array4<Real const> const& apy,
                              Array4<EBCellFlag const> const& flag,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvylo,
                              Array4<Real const> const& bvyhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apy(i,j,0) == 0.) {
                fy(i,j,0,0) = Real(0.0);
                fy(i,j,0,1) = Real(0.0);
                fy(i,j,0,2) = Real(0.0);
                fy(i,j,0,3) = Real(0.0);
            }
            else
            {
                int ihip = i + flag(i,j  ,0).isConnected( 1,0,0);
                int ihim = i - flag(i,j  ,0).isConnected(-1,0,0);
                int ilop = i + flag(i,j-1,0).isConnected( 1,0,0);
                int ilom = i - flag(i,j-1,0).isConnected(-1,0,0);
                Real whi = mlebtensor_weight(ihip-ihim);
                Real wlo = mlebtensor_weight(ilop-ilom);
                Real dudx = mlebtensor_dx_on_yface(i,j,k,0,vel,dxi,
                                                   bvylo,bvyhi,bct,dlo,dhi,
                                                   whi,wlo,ihip,ihim,ilop,ilom);
                Real dvdx = mlebtensor_dx_on_yface(i,j,k,1,vel,dxi,
                                                   bvylo,bvyhi,bct,dlo,dhi,
                                                   whi,wlo,ihip,ihim,ilop,ilom);

                Real dudy = (vel(i,j,0,0) - vel(i,j-1,0,0))*dyi;
                Real dvdy = (vel(i,j,0,1) - vel(i,j-1,0,1))*dyi;
                fy(i,j,0,0) = dudx;
                fy(i,j,0,1) = dvdx;
                fy(i,j,0,2) = dudy;
                fy(i,j,0,3) = dvdy;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel, Array4<Real const> const& apx,
                              Array4<EBCellFlag const> const& flag,
                              Array4<int const> const& ccm,
                              Array4<Real const> const& fcx,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apx(i,j,0) == Real(0.0))
            {
                fx(i,j,0,0) = Real(0.0);
                fx(i,j,0,1) = Real(0.0);
                fx(i,j,0,2) = Real(0.0);
                fx(i,j,0,3) = Real(0.0);
            }
            else
            {
                Real dudx =  (vel(i,j,0,0) - vel(i-1,j,0,0))*dxi;
                Real dvdx =  (vel(i,j,0,1) - vel(i-1,j,0,1))*dxi;
                if (apx(i,j,0) < Real(1.0)) {
                    int jj = j + (fcx(i,j,0,0) >= Real(0.0) ? 1 : -1);
                    Real fracy = (ccm(i-1,jj,0) || ccm(i,jj,0)) ? std::abs(fcx(i,j,0,0)) : Real(0.0);
                    dudx =  (Real(1.0)-fracy)*dudx
                          + fracy      *(vel(i,jj,0,0) - vel(i-1,jj,0,0))*dxi;
                    dvdx =  (Real(1.0)-fracy)*dvdx
                          + fracy      *(vel(i,jj,0,1) - vel(i-1,jj,0,1))*dxi;
                }

                int jhip = j + flag(i  ,j,0).isConnected(0, 1,0);
                int jhim = j - flag(i  ,j,0).isConnected(0,-1,0);
                int jlop = j + flag(i-1,j,0).isConnected(0, 1,0);
                int jlom = j - flag(i-1,j,0).isConnected(0,-1,0);
                Real whi = mlebtensor_weight(jhip-jhim);
                Real wlo = mlebtensor_weight(jlop-jlom);
                Real dudy = mlebtensor_dy_on_xface(i,j,k,0,vel,dyi,
                                                   whi,wlo,jhip,jhim,jlop,jlom);
                Real dvdy = mlebtensor_dy_on_xface(i,j,k,1,vel,dyi,
                                                   whi,wlo,jhip,jhim,jlop,jlom);
                fx(i,j,0,0) = dudx;
                fx(i,j,0,1) = dvdx;
                fx(i,j,0,2) = dudy;
                fx(i,j,0,3) = dvdy;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel, Array4<Real const> const& apy,
                              Array4<EBCellFlag const> const& flag,
                              Array4<int const> const& ccm,
                              Array4<Real const> const& fcy,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apy(i,j,0) == 0.) {
                fy(i,j,0,0) = Real(0.0);
                fy(i,j,0,1) = Real(0.0);
                fy(i,j,0,2) = Real(0.0);
                fy(i,j,0,3) = Real(0.0);
            }
            else
            {
                int ihip = i + flag(i,j  ,0).isConnected( 1,0,0);
                int ihim = i - flag(i,j  ,0).isConnected(-1,0,0);
                int ilop = i + flag(i,j-1,0).isConnected( 1,0,0);
                int ilom = i - flag(i,j-1,0).isConnected(-1,0,0);
                Real whi = mlebtensor_weight(ihip-ihim);
                Real wlo = mlebtensor_weight(ilop-ilom);
                Real dudx = mlebtensor_dx_on_yface(i,j,k,0,vel,dxi,
                                                   whi,wlo,ihip,ihim,ilop,ilom);
                Real dvdx = mlebtensor_dx_on_yface(i,j,k,1,vel,dxi,
                                                   whi,wlo,ihip,ihim,ilop,ilom);

                Real dudy =  (vel(i,j,0,0) - vel(i,j-1,0,0))*dyi;
                Real dvdy =  (vel(i,j,0,1) - vel(i,j-1,0,1))*dyi;
                if (apy(i,j,0) < Real(1.0)) {
                    int ii = i + (fcy(i,j,0,0) >= Real(0.0) ? 1 : -1);
                    Real fracx = (ccm(ii,j-1,0) || ccm(ii,j,0)) ? std::abs(fcy(i,j,0,0)) : Real(0.0);
                    dudy =  (Real(1.0)-fracx)*dudy
                          + fracx      *(vel(ii,j,0,0) - vel(ii,j-1,0,0))*dyi;
                    dvdy =  (Real(1.0)-fracx)*dvdy
                          + fracx      *(vel(ii,j,0,1) - vel(ii,j-1,0,1))*dyi;
                }

                fy(i,j,0,0) = dudx;
                fy(i,j,0,1) = dvdx;
                fy(i,j,0,2) = dudy;
                fy(i,j,0,3) = dvdy;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel, Array4<Real const> const& apx,
                              Array4<EBCellFlag const> const& flag,
                              Array4<int const> const& ccm,
                              Array4<Real const> const& fcx,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvxlo,
                              Array4<Real const> const& bvxhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apx(i,j,0) == Real(0.0))
            {
                fx(i,j,0,0) = Real(0.0);
                fx(i,j,0,1) = Real(0.0);
                fx(i,j,0,2) = Real(0.0);
                fx(i,j,0,3) = Real(0.0);
            }
            else
            {
                Real dudx = (vel(i,j,0,0) - vel(i-1,j,0,0))*dxi;
                Real dvdx = (vel(i,j,0,1) - vel(i-1,j,0,1))*dxi;
                if (apx(i,j,0) < Real(1.0)) {
                    int jj = j + (fcx(i,j,0,0) >= Real(0.0) ? 1 : -1);
                    Real fracy = (ccm(i-1,jj,0) || ccm(i,jj,0)) ? std::abs(fcx(i,j,0,0)) : Real(0.0);
                    dudx =  (Real(1.0)-fracy)*dudx
                          + fracy      *(vel(i,jj,0,0) - vel(i-1,jj,0,0))*dxi;
                    dvdx =  (Real(1.0)-fracy)*dvdx
                          + fracy      *(vel(i,jj,0,1) - vel(i-1,jj,0,1))*dxi;
                }

                int jhip = j + flag(i  ,j,0).isConnected(0, 1,0);
                int jhim = j - flag(i  ,j,0).isConnected(0,-1,0);
                int jlop = j + flag(i-1,j,0).isConnected(0, 1,0);
                int jlom = j - flag(i-1,j,0).isConnected(0,-1,0);
                Real whi = mlebtensor_weight(jhip-jhim);
                Real wlo = mlebtensor_weight(jlop-jlom);
                Real dudy = mlebtensor_dy_on_xface(i,j,k,0,vel,dyi,
                                                   bvxlo,bvxhi,bct,dlo,dhi,
                                                   whi,wlo,jhip,jhim,jlop,jlom);
                Real dvdy = mlebtensor_dy_on_xface(i,j,k,1,vel,dyi,
                                                   bvxlo,bvxhi,bct,dlo,dhi,
                                                   whi,wlo,jhip,jhim,jlop,jlom);
                fx(i,j,0,0) = dudx;
                fx(i,j,0,1) = dvdx;
                fx(i,j,0,2) = dudy;
                fx(i,j,0,3) = dvdy;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_vel_grads_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel, Array4<Real const> const& apy,
                              Array4<EBCellFlag const> const& flag,
                              Array4<int const> const& ccm,
                              Array4<Real const> const& fcy,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvylo,
                              Array4<Real const> const& bvyhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apy(i,j,0) == 0.) {
                fy(i,j,0,0) = Real(0.0);
                fy(i,j,0,1) = Real(0.0);
                fy(i,j,0,2) = Real(0.0);
                fy(i,j,0,3) = Real(0.0);
            }
            else
            {
                int ihip = i + flag(i,j  ,0).isConnected( 1,0,0);
                int ihim = i - flag(i,j  ,0).isConnected(-1,0,0);
                int ilop = i + flag(i,j-1,0).isConnected( 1,0,0);
                int ilom = i - flag(i,j-1,0).isConnected(-1,0,0);
                Real whi = mlebtensor_weight(ihip-ihim);
                Real wlo = mlebtensor_weight(ilop-ilom);
                Real dudx = mlebtensor_dx_on_yface(i,j,k,0,vel,dxi,
                                                   bvylo,bvyhi,bct,dlo,dhi,
                                                   whi,wlo,ihip,ihim,ilop,ilom);
                Real dvdx = mlebtensor_dx_on_yface(i,j,k,1,vel,dxi,
                                                   bvylo,bvyhi,bct,dlo,dhi,
                                                   whi,wlo,ihip,ihim,ilop,ilom);

                Real dudy = (vel(i,j,0,0) - vel(i,j-1,0,0))*dyi;
                Real dvdy = (vel(i,j,0,1) - vel(i,j-1,0,1))*dyi;
                if (apy(i,j,0) < Real(1.0)) {
                    int ii = i + (fcy(i,j,0,0) >= Real(0.0) ? 1 : -1);
                    Real fracx = (ccm(ii,j-1,0) || ccm(ii,j,0)) ? std::abs(fcy(i,j,0,0)) : Real(0.0);
                    dudy =  (Real(1.0)-fracx)*dudy
                          + fracx      *(vel(ii,j,0,0) - vel(ii,j-1,0,0))*dyi;
                    dvdy =  (Real(1.0)-fracx)*dvdy
                          + fracx      *(vel(ii,j,0,1) - vel(ii,j-1,0,1))*dyi;
                }

                fy(i,j,0,0) = dudx;
                fy(i,j,0,1) = dvdx;
                fy(i,j,0,2) = dudy;
                fy(i,j,0,3) = dvdy;
            }
        }
    }
}

}
#endif
