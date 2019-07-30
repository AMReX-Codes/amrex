#ifndef AMREX_ML_EB_TENSOR_3D_K_H_
#define AMREX_ML_EB_TENSOR_3D_K_H_

#include <AMReX_MLEBABecLap_K.H>

namespace amrex {

namespace {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real mlebtensor_weight (int d) {
        return (d==2) ? 0.5 : ((d==1) ? 1.0 : 0.0);
    }
}

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
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apx(i,j,k) == 0.0)
                {
                    fx(i,j,k,0) = 0.0;
                    fx(i,j,k,1) = 0.0;
                    fx(i,j,k,2) = 0.0;
                }
                else
                {
                    int jhip = j + flag(i  ,j,k).isConnected(0, 1,0);
                    int jhim = j - flag(i  ,j,k).isConnected(0,-1,0);
                    int jlop = j + flag(i-1,j,k).isConnected(0, 1,0);
                    int jlom = j - flag(i-1,j,k).isConnected(0,-1,0);
                    Real whi = mlebtensor_weight(jhip-jhim);
                    Real wlo = mlebtensor_weight(jlop-jlom);
                    Real dudy = (0.5*dyi) * ((vel(i  ,jhip,k,0)-vel(i  ,jhim,k,0))*whi
                                            +(vel(i-1,jlop,k,0)-vel(i-1,jlom,k,0))*wlo);
                    Real dvdy = (0.5*dyi) * ((vel(i  ,jhip,k,1)-vel(i  ,jhim,k,1))*whi
                                            +(vel(i-1,jlop,k,1)-vel(i-1,jlom,k,1))*wlo);

                    int khip = k + flag(i  ,j,k).isConnected(0,0, 1);
                    int khim = k - flag(i  ,j,k).isConnected(0,0,-1);
                    int klop = k + flag(i-1,j,k).isConnected(0,0, 1);
                    int klom = k - flag(i-1,j,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dudz = (0.5*dzi) * ((vel(i  ,j,khip,0)-vel(i  ,j,khim,0))*whi
                                            +(vel(i-1,j,klop,0)-vel(i-1,j,klom,0))*wlo);
                    Real dwdz = (0.5*dzi) * ((vel(i  ,j,khip,2)-vel(i  ,j,khim,2))*whi
                                            +(vel(i-1,j,klop,2)-vel(i-1,j,klom,2))*wlo);

                    Real divu = dvdy + dwdz;
                    Real xif = kapx(i,j,k);
                    Real mun = 0.75*(etax(i,j,k,0)-xif);  // restore the original eta
                    Real mut =       etax(i,j,k,1);
                    fx(i,j,k,0) = -mun*(-twoThirds*divu) - xif*divu;
                    fx(i,j,k,1) = -mut*dudy;
                    fx(i,j,k,2) = -mut*dudz;
                }
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
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apy(i,j,k) == 0.0)
                {
                    fy(i,j,k,0) = 0.0;
                    fy(i,j,k,1) = 0.0;
                    fy(i,j,k,2) = 0.0;
                }
                else
                {
                    int ihip = i + flag(i,j  ,k).isConnected( 1,0,0);
                    int ihim = i - flag(i,j  ,k).isConnected(-1,0,0);
                    int ilop = i + flag(i,j-1,k).isConnected( 1,0,0);
                    int ilom = i - flag(i,j-1,k).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);
                    Real dudx = (0.5*dxi) * ((vel(ihip,j  ,k,0)-vel(ihim,j  ,k,0))*whi
                                            +(vel(ilop,j-1,k,0)-vel(ilom,j-1,k,0))*wlo);
                    Real dvdx = (0.5*dxi) * ((vel(ihip,j  ,k,1)-vel(ihim,j  ,k,1))*whi
                                            +(vel(ilop,j-1,k,1)-vel(ilom,j-1,k,1))*wlo);

                    int khip = k + flag(i,j  ,k).isConnected(0,0, 1);
                    int khim = k - flag(i,j  ,k).isConnected(0,0,-1);
                    int klop = k + flag(i,j-1,k).isConnected(0,0, 1);
                    int klom = k - flag(i,j-1,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dvdz = (0.5*dzi) * ((vel(i,j  ,khip,1)-vel(i,j  ,khim,1))*whi
                                            +(vel(i,j-1,klop,1)-vel(i,j-1,klom,1))*wlo);
                    Real dwdz = (0.5*dzi) * ((vel(i,j  ,khip,2)-vel(i,j  ,khim,2))*whi
                                            +(vel(i,j-1,klop,2)-vel(i,j-1,klom,2))*wlo);

                    Real divu = dudx + dwdz;
                    Real xif = kapy(i,j,k);
                    Real mun = 0.75*(etay(i,j,k,1)-xif);  // restore the original eta
                    Real mut =       etay(i,j,k,0);
                    fy(i,j,k,0) = -mut*dvdx;
                    fy(i,j,k,1) = -mun*(-twoThirds*divu) - xif*divu;
                    fy(i,j,k,2) = -mut*dvdz;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms_fz (Box const& box, Array4<Real> const& fz,
                                Array4<Real const> const& vel,
                                Array4<Real const> const& etaz,
                                Array4<Real const> const& kapz,
                                Array4<Real const> const& apz,
                                Array4<EBCellFlag const> const& flag,
                                GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (apz(i,j,k) == 0.0)
                {
                    fz(i,j,k,0) = 0.0;
                    fz(i,j,k,1) = 0.0;
                    fz(i,j,k,2) = 0.0;
                }
                else
                {
                    int ihip = i + flag(i,j,k  ).isConnected( 1,0,0);
                    int ihim = i - flag(i,j,k  ).isConnected(-1,0,0);
                    int ilop = i + flag(i,j,k-1).isConnected( 1,0,0);
                    int ilom = i - flag(i,j,k-1).isConnected(-1,0,0);
                    Real whi = mlebtensor_weight(ihip-ihim);
                    Real wlo = mlebtensor_weight(ilop-ilom);

                    Real dudx = (0.5*dxi) * ((vel(ihip,j,k  ,0)-vel(ihim,j,k  ,0))*whi
                                            +(vel(ilop,j,k-1,0)-vel(ilom,j,k-1,0))*wlo);
                    Real dwdx = (0.5*dxi) * ((vel(ihip,j,k  ,2)-vel(ihim,j,k  ,2))*whi
                                            +(vel(ilop,j,k-1,2)-vel(ilom,j,k-1,2))*wlo);

                    int jhip = j + flag(i,j,k  ).isConnected(0, 1,0);
                    int jhim = j - flag(i,j,k  ).isConnected(0,-1,0);
                    int jlop = j + flag(i,j,k-1).isConnected(0, 1,0);
                    int jlom = j - flag(i,j,k-1).isConnected(0,-1,0);
                    whi = mlebtensor_weight(jhip-jhim);
                    wlo = mlebtensor_weight(jlop-jlom);
                    Real dvdy = (0.5*dyi) * ((vel(i,jhip,k  ,1)-vel(i,jhim,k  ,1))*whi
                                            +(vel(i,jlop,k-1,1)-vel(i,jlom,k-1,1))*wlo);
                    Real dwdy = (0.5*dyi) * ((vel(i,jhip,k  ,2)-vel(i,jhim,k  ,2))*whi
                                            +(vel(i,jlop,k-1,2)-vel(i,jlom,k-1,2))*wlo);

                    Real divu = dudx + dvdy;
                    Real xif = kapz(i,j,k);
                    Real mun = 0.75*(etaz(i,j,k,2)-xif);  // restore the original eta
                    Real mut =       etaz(i,j,k,0);

                    fz(i,j,k,0) = -mut*dwdx;
                    fz(i,j,k,1) = -mut*dwdy;
                    fz(i,j,k,2) = -mun*(-twoThirds*divu) - xif*divu;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebtensor_cross_terms (Box const& box, Array4<Real> const& Ax,
                             Array4<Real const> const& fx,
                             Array4<Real const> const& fy,
                             Array4<Real const> const& fz,
                             Array4<Real const> const& vel,
//                             Array4<Real const> const& velb,
                             Array4<Real const> const& etab,
                             Array4<Real const> const& kapb,
                             Array4<int const> const& ccm,
                             Array4<EBCellFlag const> const& flag,
                             Array4<Real const> const& vol,
                             Array4<Real const> const& apx,
                             Array4<Real const> const& apy,
                             Array4<Real const> const& apz,
                             Array4<Real const> const& fcx,
                             Array4<Real const> const& fcy,
                             Array4<Real const> const& fcz,
                             Array4<Real const> const& bc,
//                             int is_dirichlet, int is_inhomog,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                             Real bscalar) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (flag(i,j,k).isRegular())
                {
                    Ax(i,j,k,0) += bscalar*(dxi*(fx(i+1,j  ,k  ,0) - fx(i,j,k,0))
                                          + dyi*(fy(i  ,j+1,k  ,0) - fy(i,j,k,0))
                                          + dzi*(fz(i  ,j  ,k+1,0) - fz(i,j,k,0)));
                    Ax(i,j,k,1) += bscalar*(dxi*(fx(i+1,j  ,k  ,1) - fx(i,j,k,1))
                                          + dyi*(fy(i  ,j+1,k  ,1) - fy(i,j,k,1))
                                          + dzi*(fz(i  ,j  ,k+1,1) - fz(i,j,k,1)));
                    Ax(i,j,k,2) += bscalar*(dxi*(fx(i+1,j  ,k  ,2) - fx(i,j,k,2))
                                          + dyi*(fy(i  ,j+1,k  ,2) - fy(i,j,k,2))
                                          + dzi*(fz(i  ,j  ,k+1,2) - fz(i,j,k,2)));
                }
                else if (flag(i,j,k).isSingleValued())
                {
                    Real fxm_0 = fx(i,j,k,0);
                    Real fxm_1 = fx(i,j,k,1);
                    Real fxm_2 = fx(i,j,k,2);
                    if (apx(i,j,k) > 0.0 and apx(i,j,k) < 1.0) {
                        int jj = j + (fcx(i,j,k,0) >= 0.0 ? 1 : -1);
                        int kk = k + (fcx(i,j,k,1) >= 0.0 ? 1 : -1);
                        Real fracy = (ccm(i-1,jj,k) or ccm(i,jj,k)) ? std::abs(fcx(i,j,k,0)) : 0.0;
                        Real fracz = (ccm(i-1,j,kk) or ccm(i,j,kk)) ? std::abs(fcx(i,j,k,1)) : 0.0;
                        fxm_0 = (1.0-fracy)*(1.0-fracz) *fxm_0
                            + fracy*(1.0-fracz) * fx(i,jj,k ,0)
                            + fracz*(1.0-fracy) * fx(i,j ,kk,0)
                            + fracy*     fracz  * fx(i,jj,kk,0);
                        fxm_1 = (1.0-fracy)*(1.0-fracz) *fxm_1
                            + fracy*(1.0-fracz) * fx(i,jj,k ,1)
                            + fracz*(1.0-fracy) * fx(i,j ,kk,1)
                            + fracy*     fracz  * fx(i,jj,kk,1);
                        fxm_2 = (1.0-fracy)*(1.0-fracz) *fxm_2
                            + fracy*(1.0-fracz) * fx(i,jj,k ,2)
                            + fracz*(1.0-fracy) * fx(i,j ,kk,2)
                            + fracy*     fracz  * fx(i,jj,kk,2);
                    }

                    Real fxp_0 = fx(i+1,j,k,0);
                    Real fxp_1 = fx(i+1,j,k,1);
                    Real fxp_2 = fx(i+1,j,k,2);
                    if (apx(i+1,j,k) > 0.0 and apx(i+1,j,k) < 1.0) {
                        int jj = j + (fcx(i+1,j,k,0) >= 0.0 ? 1 : -1);
                        int kk = k + (fcx(i+1,j,k,1) >= 0.0 ? 1 : -1);
                        Real fracy = (ccm(i,jj,k) or ccm(i+1,jj,k)) ? std::abs(fcx(i+1,j,k,0)) : 0.0;
                        Real fracz = (ccm(i,j,kk) or ccm(i+1,j,kk)) ? std::abs(fcx(i+1,j,k,1)) : 0.0;
                        fxp_0 = (1.0-fracy)*(1.0-fracz) *  fxp_0
                            + fracy*(1.0-fracz) * fx(i+1,jj,k ,0)
                            + fracz*(1.0-fracy) * fx(i+1,j ,kk,0)
                            + fracy*     fracz  * fx(i+1,jj,kk,0);
                        fxp_1 = (1.0-fracy)*(1.0-fracz) *fxp_1
                            + fracy*(1.0-fracz) * fx(i+1,jj,k ,1)
                            + fracz*(1.0-fracy) * fx(i+1,j ,kk,1)
                            + fracy*     fracz  * fx(i+1,jj,kk,1);
                        fxp_2 = (1.0-fracy)*(1.0-fracz) *fxp_2
                            + fracy*(1.0-fracz) * fx(i+1,jj,k ,2)
                            + fracz*(1.0-fracy) * fx(i+1,j ,kk,2)
                            + fracy*     fracz  * fx(i+1,jj,kk,2);
                    }

                    Real fym_0 = fy(i,j,k,0);
                    Real fym_1 = fy(i,j,k,1);
                    Real fym_2 = fy(i,j,k,2);
                    if (apy(i,j,k) > 0.0 and apy(i,j,k) < 1.0) {
                        int ii = i + (fcy(i,j,k,0) >= 0.0 ? 1 : -1);
                        int kk = k + (fcy(i,j,k,1) >= 0.0 ? 1 : -1);
                        Real fracx = (ccm(ii,j-1,k) or ccm(ii,j,k)) ? std::abs(fcy(i,j,k,0)) : 0.0;
                        Real fracz = (ccm(i,j-1,kk) or ccm(i,j,kk)) ? std::abs(fcy(i,j,k,1)) : 0.0;
                        fym_0 = (1.0-fracx)*(1.0-fracz) *fym_0
                            + fracx*(1.0-fracz) * fy(ii,j,k ,0)
                            + fracz*(1.0-fracx) * fy(i ,j,kk,0)
                            + fracx*     fracz  * fy(ii,j,kk,0);
                        fym_1 = (1.0-fracx)*(1.0-fracz) *fym_1
                            + fracx*(1.0-fracz) * fy(ii,j,k ,1)
                            + fracz*(1.0-fracx) * fy(i ,j,kk,1)
                            + fracx*     fracz  * fy(ii,j,kk,1);
                        fym_2 = (1.0-fracx)*(1.0-fracz) *fym_2
                            + fracx*(1.0-fracz) * fy(ii,j,k ,2)
                            + fracz*(1.0-fracx) * fy(i ,j,kk,2)
                            + fracx*     fracz  * fy(ii,j,kk,2);
                    }

                    Real fyp_0 = fy(i,j+1,k,0);
                    Real fyp_1 = fy(i,j+1,k,1);
                    Real fyp_2 = fy(i,j+1,k,2);
                    if (apy(i,j+1,k) > 0.0 and apy(i,j+1,k) < 1.0) {
                        int ii = i + (fcy(i,j+1,k,0) >= 0.0 ? 1 : -1);
                        int kk = k + (fcy(i,j+1,k,1) >= 0.0 ? 1 : -1);
                        Real fracx = (ccm(ii,j,k) or ccm(ii,j+1,k)) ? std::abs(fcy(i,j+1,k,0)) : 0.0;
                        Real fracz = (ccm(i,j,kk) or ccm(i,j+1,kk)) ? std::abs(fcy(i,j+1,k,1)) : 0.0;
                        fyp_0 = (1.0-fracx)*(1.0-fracz) *  fyp_0
                            + fracx*(1.0-fracz) * fy(ii,j+1,k ,0)
                            + fracz*(1.0-fracx) * fy(i ,j+1,kk,0)
                            + fracx*     fracz  * fy(ii,j+1,kk,0);
                        fyp_1 = (1.0-fracx)*(1.0-fracz) *fyp_1
                            + fracx*(1.0-fracz) * fy(ii,j+1,k ,1)
                            + fracz*(1.0-fracx) * fy(i ,j+1,kk,1)
                            + fracx*     fracz  * fy(ii,j+1,kk,1);
                        fyp_2 = (1.0-fracx)*(1.0-fracz) *fyp_2
                            + fracx*(1.0-fracz) * fy(ii,j+1,k ,2)
                            + fracz*(1.0-fracx) * fy(i ,j+1,kk,2)
                            + fracx*     fracz  * fy(ii,j+1,kk,2);
                    }

                    Real fzm_0 = fz(i,j,k,0);
                    Real fzm_1 = fz(i,j,k,1);
                    Real fzm_2 = fz(i,j,k,2);
                    if (apz(i,j,k) > 0.0 and apz(i,j,k) < 1.0) {
                        int ii = i + (fcz(i,j,k,0) >= 0.0 ? 1 : -1);
                        int jj = j + (fcz(i,j,k,1) >= 0.0 ? 1 : -1);
                        Real fracx = (ccm(ii,j,k-1) or ccm(ii,j,k)) ? std::abs(fcz(i,j,k,0)) : 0.0;
                        Real fracy = (ccm(i,jj,k-1) or ccm(i,jj,k)) ? std::abs(fcz(i,j,k,1)) : 0.0;
                        fzm_0 = (1.0-fracx)*(1.0-fracy) *fzm_0
                            + fracx*(1.0-fracy) * fz(ii,j ,k,0)
                            + fracy*(1.0-fracx) * fz(i ,jj,k,0)
                            + fracx*     fracy  * fz(ii,jj,k,0);
                        fzm_1 = (1.0-fracx)*(1.0-fracy) *fzm_1
                            + fracx*(1.0-fracy) * fz(ii,j ,k,1)
                            + fracy*(1.0-fracx) * fz(i ,jj,k,1)
                            + fracx*     fracy  * fz(ii,jj,k,1);
                        fzm_2 = (1.0-fracx)*(1.0-fracy) *fzm_2
                            + fracx*(1.0-fracy) * fz(ii,j ,k,2)
                            + fracy*(1.0-fracx) * fz(i ,jj,k,2)
                            + fracx*     fracy  * fz(ii,jj,k,2);
                    }

                    Real fzp_0 = fz(i,j,k+1,0);
                    Real fzp_1 = fz(i,j,k+1,1);
                    Real fzp_2 = fz(i,j,k+1,2);
                    if (apz(i,j,k+1) > 0.0 and apz(i,j,k+1) < 1.0) {
                        int ii = i + (fcz(i,j,k+1,0) >= 0.0 ? 1 : -1);
                        int jj = j + (fcz(i,j,k+1,1) >= 0.0 ? 1 : -1);
                        Real fracx = (ccm(ii,j,k) or ccm(ii,j,k)) ? std::abs(fcz(i,j,k+1,0)) : 0.0;
                        Real fracy = (ccm(i,jj,k) or ccm(i,jj,k)) ? std::abs(fcz(i,j,k+1,1)) : 0.0;
                        fzp_0 = (1.0-fracx)*(1.0-fracy) *  fzp_0
                            + fracx*(1.0-fracy) * fz(ii,j ,k+1,0)
                            + fracy*(1.0-fracx) * fz(i ,jj,k+1,0)
                            + fracx*     fracy  * fz(ii,jj,k+1,0);
                        fzp_1 = (1.0-fracx)*(1.0-fracy) *fzp_1
                            + fracx*(1.0-fracy) * fz(ii,j ,k+1,1)
                            + fracy*(1.0-fracx) * fz(i ,jj,k+1,1)
                            + fracx*     fracy  * fz(ii,jj,k+1,1);
                        fzp_2 = (1.0-fracx)*(1.0-fracy) *fzp_2
                            + fracx*(1.0-fracy) * fz(ii,j ,k+1,2)
                            + fracy*(1.0-fracx) * fz(i ,jj,k+1,2)
                            + fracx*     fracy  * fz(ii,jj,k+1,2);
                    }

                    Real kappa = vol(i,j,k);
                    Real feb_0 = 0.0;
                    Real feb_1 = 0.0;
                    Real feb_2 = 0.0;
                    // if (is_dirichlet) 
                    {
                        Real dapx = apx(i+1,j,k)-apx(i,j,k);
                        Real dapy = apy(i,j+1,k)-apy(i,j,k);
                        Real dapz = apz(i,j,k+1)-apz(i,j,k);
                        Real anorminv = 1.0/std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
                        Real anrmx = -dapx * anorminv;
                        Real anrmy = -dapy * anorminv;
                        Real anrmz = -dapz * anorminv;

                        Real velb_0 = 0, velb_1 = 0, velb_2 = 0;
#if 0
                        if (is_inhomog) {
                            velb_0 = velb(i,j,k,0);
                            velb_1 = velb(i,j,k,1);
                            velb_2 = velb(i,j,k,2);
                        }
#endif

                        Real dx_eb = amrex::max(0.3, (kappa*kappa-0.25)/(2.*kappa));
                        Real dg = dx_eb / amrex::max(std::abs(anrmx),
                                                     std::abs(anrmy),
                                                     std::abs(anrmz));
                        Real dginv = 1.0/dg;
                        Real gx = bc(i,j,k,0) - dg*anrmx;
                        Real gy = bc(i,j,k,1) - dg*anrmy;
                        Real gz = bc(i,j,k,2) - dg*anrmz;
                        int isx = (anrmx > 0.0) ? 1 : -1;
                        int isy = (anrmy > 0.0) ? 1 : -1;
                        int isz = (anrmz > 0.0) ? 1 : -1;
                        int ii = i - isx;
                        int jj = j - isy;
                        int kk = k - isz;

                        gx *= isx;
                        gy *= isy;
                        gz *= isz;
                        Real gxy = gx*gy;
                        Real gxz = gx*gz;
                        Real gyz = gy*gz;
                        Real gxyz = gx*gy*gz;
                        Real oneggg = 1.0+gx+gy+gz+gxy+gxz+gyz+gxyz;

                        Real velg = oneggg             * vel(i ,j ,k ,0)
                            + (-gz - gxz - gyz - gxyz) * vel(i ,j ,kk,0)
                            + (-gy - gxy - gyz - gxyz) * vel(i ,jj,k ,0)
                            + (gyz + gxyz)             * vel(i ,jj,kk,0)
                            + (-gx - gxy - gxz - gxyz) * vel(ii,j ,k ,0)
                            + (gxz + gxyz)             * vel(ii,j ,kk,0)
                            + (gxy + gxyz)             * vel(ii,jj,k ,0)
                            + (-gxyz)                  * vel(ii,jj,kk,0);
                        Real dudn = (velb_0-velg) * dginv;

                        velg      = oneggg             * vel(i ,j ,k ,1)
                            + (-gz - gxz - gyz - gxyz) * vel(i ,j ,kk,1)
                            + (-gy - gxy - gyz - gxyz) * vel(i ,jj,k ,1)
                            + (gyz + gxyz)             * vel(i ,jj,kk,1)
                            + (-gx - gxy - gxz - gxyz) * vel(ii,j ,k ,1)
                            + (gxz + gxyz)             * vel(ii,j ,kk,1)
                            + (gxy + gxyz)             * vel(ii,jj,k ,1)
                            + (-gxyz)                  * vel(ii,jj,kk,1);
                        Real dvdn = (velb_1-velg) * dginv;

                        velg      = oneggg             * vel(i ,j ,k ,2)
                            + (-gz - gxz - gyz - gxyz) * vel(i ,j ,kk,2)
                            + (-gy - gxy - gyz - gxyz) * vel(i ,jj,k ,2)
                            + (gyz + gxyz)             * vel(i ,jj,kk,2)
                            + (-gx - gxy - gxz - gxyz) * vel(ii,j ,k ,2)
                            + (gxz + gxyz)             * vel(ii,j ,kk,2)
                            + (gxy + gxyz)             * vel(ii,jj,k ,2)
                            + (-gxyz)                  * vel(ii,jj,kk,2);
                        Real dwdn = (velb_2-velg) * dginv;

                        // transverse derivatives are zero on the wall
                        Real dudx = dudn * anrmx;
                        Real dudy = dudn * anrmy;
                        Real dudz = dudn * anrmz;
                        Real dvdx = dvdn * anrmx;
                        Real dvdy = dvdn * anrmy;
                        Real dvdz = dvdn * anrmz;
                        Real dwdx = dwdn * anrmx;
                        Real dwdy = dwdn * anrmy;
                        Real dwdz = dwdn * anrmz;
                        Real divu = dudx + dvdy + dwdz;
                        Real xi = kapb(i,j,k);
                        Real mu = etab(i,j,k);
                        Real tautmp = (xi-(2./3.)*mu)*divu;
                        // Note that mu*(grad v) has been included already in MLEBABecLap
                        Real tauxx = mu*dudx + tautmp;
                        Real tauyx = mu*dvdx;
                        Real tauzx = mu*dwdx;
                        Real tauyy = mu*dvdy + tautmp;
                        Real tauxy = mu*dudy;
                        Real tauzy = mu*dwdy;
                        Real tauzz = mu*dwdz + tautmp;
                        Real tauxz = mu*dudz;
                        Real tauyz = mu*dvdz;
                        // assuming dxi == dyi == dzi
                        feb_0 = dxi*(dapx*tauxx + dapy*tauyx + dapz*tauzx);
                        feb_1 = dxi*(dapx*tauxy + dapy*tauyy + dapz*tauzy);
                        feb_2 = dxi*(dapx*tauxz + dapy*tauyz + dapz*tauzz);
                    }

                    Real volinv = bscalar / kappa;
                    Ax(i,j,k,0) += volinv * (dxi*(apx(i+1,j,k)*fxp_0-apx(i,j,k)*fxm_0)
                                            +dyi*(apy(i,j+1,k)*fyp_0-apy(i,j,k)*fym_0)
                                            +dzi*(apz(i,j,k+1)*fzp_0-apz(i,j,k)*fzm_0)
                                            +dxi*                               feb_0);
                    Ax(i,j,k,1) += volinv * (dxi*(apx(i+1,j,k)*fxp_1-apx(i,j,k)*fxm_1)
                                            +dyi*(apy(i,j+1,k)*fyp_1-apy(i,j,k)*fym_1)
                                            +dzi*(apz(i,j,k+1)*fzp_1-apz(i,j,k)*fzm_1)
                                            +dxi*                               feb_1);
                    Ax(i,j,k,2) += volinv * (dxi*(apx(i+1,j,k)*fxp_2-apx(i,j,k)*fxm_2)
                                            +dyi*(apy(i,j+1,k)*fyp_2-apy(i,j,k)*fym_2)
                                            +dzi*(apz(i,j,k+1)*fzp_2-apz(i,j,k)*fzm_2)
                                            +dxi*                               feb_2);
                }
            }
        }
    }
}

}

#endif
