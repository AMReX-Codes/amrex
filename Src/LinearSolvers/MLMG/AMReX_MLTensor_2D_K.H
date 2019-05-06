#ifndef AMREX_MLTENSOR_2D_K_H_
#define AMREX_MLTENSOR_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mltensor_cross_terms_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etax,
                              Array4<Real const> const& kapx,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudx = (vel(i,j,0,0)-vel(i-1,j,0,0))*dxi;
            Real dudy = (vel(i,j+1,0,0)+vel(i-1,j+1,0,0)-vel(i,j-1,0,0)-vel(i-1,j-1,0,0))*(0.25*dyi);
            Real dvdy = (vel(i,j+1,0,1)+vel(i-1,j+1,0,1)-vel(i,j-1,0,1)-vel(i-1,j-1,0,1))*(0.25*dyi);
            Real divu = dudx + dvdy;
            Real muf = etax(i,j,0);
            Real xif = kapx(i,j,0);
            fx(i,j,0,0) = -muf*(dudx-twoThirds*divu) - xif*divu;
            fx(i,j,0,1) = -muf*dudy;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mltensor_cross_terms_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etay,
                              Array4<Real const> const& kapy,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dvdy = (vel(i,j,0,1)-vel(i,j-1,0,1))*dyi;
            Real dudx = (vel(i+1,j,0,0)+vel(i+1,j-1,0,0)-vel(i-1,j,0,0)-vel(i-1,j-1,0,0))*(0.25*dxi);
            Real dvdx = (vel(i+1,j,0,1)+vel(i+1,j-1,0,1)-vel(i-1,j,0,1)-vel(i-1,j-1,0,1))*(0.25*dxi);
            Real divu = dudx + dvdy;
            Real muf = etay(i,j,0);
            Real xif = kapy(i,j,0);
            fy(i,j,0,0) = -muf*dvdx;
            fy(i,j,0,1) = -muf*(dvdy-twoThirds*divu) - xif*divu;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mltensor_cross_terms (Box const& box, Array4<Real> const& Ax,
                           Array4<Real const> const& fx,
                           Array4<Real const> const& fy,
                           GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Ax(i,j,0,0) += dxi*(fx(i+1,j  ,0,0) - fx(i,j,0,0))
                +          dyi*(fy(i  ,j+1,0,0) - fy(i,j,0,0));
            Ax(i,j,0,1) += dxi*(fx(i+1,j  ,0,1) - fx(i,j,0,1))
                +          dyi*(fy(i  ,j+1,0,1) - fy(i,j,0,1));
        }
    }
}

}

#endif
