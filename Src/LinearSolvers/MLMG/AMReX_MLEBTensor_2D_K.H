#ifndef AMREX_ML_EB_TENSOR_2D_K_H_
#define AMREX_ML_EB_TENSOR_2D_K_H_

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
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apx(i,j,0) == 0.0)
            {
                fx(i,j,0,0) = 0.0;
                fx(i,j,0,1) = 0.0;
            }
            else
            {
                int jhip = j + flag(i  ,j,0).isConnected(0, 1,0);
                int jhim = j - flag(i  ,j,0).isConnected(0,-1,0);
                int jlop = j + flag(i-1,j,0).isConnected(0, 1,0);
                int jlom = j - flag(i-1,j,0).isConnected(0,-1,0);
                Real whi = mlebtensor_weight(jhip-jhim);
                Real wlo = mlebtensor_weight(jlop-jlom);
                Real dudy = (0.5*dyi) * ((vel(i  ,jhip,0,0)-vel(i  ,jhim,0,0))*whi
                                        +(vel(i-1,jlop,0,0)-vel(i-1,jlom,0,0))*wlo);
                Real dvdy = (0.5*dyi) * ((vel(i  ,jhip,0,1)-vel(i  ,jhim,0,1))*whi
                                        +(vel(i-1,jlop,0,1)-vel(i-1,jlom,0,1))*wlo);

                Real divu = dvdy;
                Real xif = kapx(i,j,0);
                Real mun = 0.75*(etax(i,j,0,0)-xif);  // restore the original eta
                Real mut =       etax(i,j,0,1);
                fx(i,j,0,0) = -mun*(-twoThirds*divu) - xif*divu;
                fx(i,j,0,1) = -mut*dudy;
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
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if (apy(i,j,0) == 0.0)
            {
                fy(i,j,0,0) = 0.0;
                fy(i,j,0,1) = 0.0;
            }
            else
            {
                int ihip = i + flag(i,j  ,0).isConnected( 1,0,0);
                int ihim = i - flag(i,j  ,0).isConnected(-1,0,0);
                int ilop = i + flag(i,j-1,0).isConnected( 1,0,0);
                int ilom = i - flag(i,j-1,0).isConnected(-1,0,0);
                Real whi = mlebtensor_weight(ihip-ihim);
                Real wlo = mlebtensor_weight(ilop-ilom);
                Real dudx = (0.5*dxi) * ((vel(ihip,j  ,0,0)-vel(ihim,j  ,0,0))*whi
                                        +(vel(ilop,j-1,0,0)-vel(ilom,j-1,0,0))*wlo);
                Real dvdx = (0.5*dxi) * ((vel(ihip,j  ,0,1)-vel(ihim,j  ,0,1))*whi
                                        +(vel(ilop,j-1,0,1)-vel(ilom,j-1,0,1))*wlo);

                Real divu = dudx;
                Real xif = kapy(i,j,0);
                Real mun = 0.75*(etay(i,j,0,1)-xif);  // restore the original eta
                Real mut =       etay(i,j,0,0);
                fy(i,j,0,0) = -mut*dvdx;
                fy(i,j,0,1) = -mun*(-twoThirds*divu) - xif*divu;
            }
        }
    }
}

}

#endif
