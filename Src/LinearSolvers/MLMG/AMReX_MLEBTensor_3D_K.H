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

AMREX_GPU_HOST_DEVICE AMREX_INLINE
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
                    fx(i,j,k,1) = -mut*(dudy);
                    fx(i,j,k,2) = -mut*(dudz);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
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
                                            +(vel(ilop,j-1,k,1)-vel(ilom,j-1,k,1))*whi);

                    int khip = k + flag(i,j  ,k).isConnected(0,0, 1);
                    int khim = k - flag(i,j  ,k).isConnected(0,0,-1);
                    int klop = k + flag(i,j-1,k).isConnected(0,0, 1);
                    int klom = k - flag(i,j-1,k).isConnected(0,0,-1);
                    whi = mlebtensor_weight(khip-khim);
                    wlo = mlebtensor_weight(klop-klom);
                    Real dvdz = (0.5*dzi) * ((vel(i,j  ,khip,1)-vel(i,j  ,khim,1))*whi
                                            +(vel(i,j-1,klop,1)-vel(i,j-1,klom,1))*wlo);
                    Real dwdz = (0.5*dzi) * ((vel(i,j  ,khip,2)-vel(i,j  ,khip,2))*whi
                                            +(vel(i,j-1,klop,2)-vel(i,j-1,klom,2))*wlo);

                    Real divu = dudx + dwdz;
                    Real xif = kapy(i,j,k);
                    Real mun = 0.75*(etay(i,j,k,1)-xif);  // restore the original eta
                    Real mut =       etay(i,j,k,0);
                    fy(i,j,k,0) = -mut*(dvdx);
                    fy(i,j,k,1) = -mun*(-twoThirds*divu) - xif*divu;
                    fy(i,j,k,2) = -mut*(dvdz);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
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
                    fz(i,j,k,0) = -mut*(dwdx);
                    fz(i,j,k,1) = -mut*(dwdy);
                    fz(i,j,k,2) = -mun*(-twoThirds*divu) - xif*divu;
                }
            }
        }
    }
}

}

#endif
