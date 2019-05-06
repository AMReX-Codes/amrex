#ifndef AMREX_MLTENSOR_3D_K_H_
#define AMREX_MLTENSOR_3D_K_H_

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
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudx = (vel(i,j,k,0)-vel(i-1,j,k,0))*dxi;
                Real dudy = (vel(i,j+1,k,0)+vel(i-1,j+1,k,0)-vel(i,j-1,k,0)-vel(i-1,j-1,k,0))*(0.25*dyi);
                Real dvdy = (vel(i,j+1,k,1)+vel(i-1,j+1,k,1)-vel(i,j-1,k,1)-vel(i-1,j-1,k,1))*(0.25*dyi);
                Real dudz = (vel(i,j,k+1,0)+vel(i-1,j,k+1,0)-vel(i,j,k-1,0)-vel(i-1,j,k-1,0))*(0.25*dzi);
                Real dwdz = (vel(i,j,k+1,2)+vel(i-1,j,k+1,2)-vel(i,j,k-1,2)-vel(i-1,j,k-1,2))*(0.25*dzi);
                Real divu = dudx + dvdy + dwdz;
                Real muf = etax(i,j,k);
                Real xif = kapx(i,j,k);
                fx(i,j,k,0) = -muf*(dudx-twoThirds*divu) - xif*divu;
                fx(i,j,k,1) = -muf*(dudy);
                fx(i,j,k,2) = -muf*(dudz);
            }
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
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dvdy = (vel(i,j,k,1)-vel(i,j-1,k,1))*dyi;
                Real dudx = (vel(i+1,j,k,0)+vel(i+1,j-1,k,0)-vel(i-1,j,k,0)-vel(i-1,j-1,k,0))*(0.25*dxi);
                Real dvdx = (vel(i+1,j,k,1)+vel(i+1,j-1,k,1)-vel(i-1,j,k,1)-vel(i-1,j-1,k,1))*(0.25*dxi);
                Real dvdz = (vel(i,j,k+1,1)+vel(i,j-1,k+1,1)-vel(i,j,k-1,1)-vel(i,j-1,k-1,1))*(0.25*dzi);
                Real dwdz = (vel(i,j,k+1,2)+vel(i,j-1,k+1,2)-vel(i,j,k-1,2)-vel(i,j-1,k-1,2))*(0.25*dzi);
                Real divu = dudx + dvdy + dwdz;
                Real muf = etay(i,j,k);
                Real xif = kapy(i,j,k);
                fy(i,j,k,0) = -muf*(dvdx);
                fy(i,j,k,1) = -muf*(dvdy-twoThirds*divu) - xif*divu;
                fy(i,j,k,2) = -muf*(dvdz);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mltensor_cross_terms_fz (Box const& box, Array4<Real> const& fz,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etaz,
                              Array4<Real const> const& kapz,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dwdz = (vel(i,j,k,2)-vel(i,j,k-1,2))*dzi;
                Real dudx = (vel(i+1,j,k,0)+vel(i+1,j,k-1,0)-vel(i-1,j,k,0)-vel(i-1,j,k-1,0))*(0.25*dxi);
                Real dwdx = (vel(i+1,j,k,2)+vel(i+1,j,k-1,2)-vel(i-1,j,k,2)-vel(i-1,j,k-1,2))*(0.25*dxi);
                Real dvdy = (vel(i,j+1,k,1)+vel(i,j+1,k-1,1)-vel(i,j-1,k,1)-vel(i,j-1,k-1,1))*(0.25*dyi);
                Real dwdy = (vel(i,j+1,k,2)+vel(i,j+1,k-1,2)-vel(i,j-1,k,2)-vel(i,j-1,k-1,2))*(0.25*dyi);
                Real divu = dudx + dvdy + dwdz;
                Real muf = etaz(i,j,k);
                Real xif = kapz(i,j,k);
                fz(i,j,k,0) = -muf*(dwdx);
                fz(i,j,k,1) = -muf*(dwdy);
                fz(i,j,k,2) = -muf*(dwdz-twoThirds*divu) - xif*divu;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mltensor_cross_terms (Box const& box, Array4<Real> const& Ax,
                           Array4<Real const> const& fx,
                           Array4<Real const> const& fy,
                           Array4<Real const> const& fz,
                           GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
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
                Ax(i,j,k,0) += dxi*(fx(i+1,j  ,k  ,0) - fx(i,j,k,0))
                    +          dyi*(fy(i  ,j+1,k  ,0) - fy(i,j,k,0))
                    +          dzi*(fz(i  ,j  ,k+1,0) - fz(i,j,k,0));
                Ax(i,j,k,1) += dxi*(fx(i+1,j  ,k  ,1) - fx(i,j,k,1))
                    +          dyi*(fy(i  ,j+1,k  ,1) - fy(i,j,k,1))
                    +          dzi*(fz(i  ,j  ,k+1,1) - fz(i,j,k,1));
                Ax(i,j,k,2) += dxi*(fx(i+1,j  ,k  ,2) - fx(i,j,k,2))
                    +          dyi*(fy(i  ,j+1,k  ,2) - fy(i,j,k,2))
                    +          dzi*(fz(i  ,j  ,k+1,2) - fz(i,j,k,2));
            }
        }
    }
}

}

#endif
