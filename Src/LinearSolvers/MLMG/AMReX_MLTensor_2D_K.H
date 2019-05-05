#ifndef AMREX_MLTENSOR_2D_K_H_
#define AMREX_MLTENSOR_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mltensor_adotx_add_extra (Box const& box, Array4<Real> const& Ax,
                               Array4<Real const> const& vel,
                               Array4<Real const> const& gradeta,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
    const Real dxi = 0.5*dxinv[0];
    const Real dyi = 0.5*dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudx = dxi*(vel(i+1,j  ,0,0) - vel(i-1,j  ,0,0));
            Real dudy = dyi*(vel(i  ,j+1,0,0) - vel(i  ,j-1,0,0));
            Real dvdx = dxi*(vel(i+1,j  ,0,1) - vel(i-1,j  ,0,1));
            Real dvdy = dyi*(vel(i  ,j+1,0,1) - vel(i  ,j-1,0,1));
            Real div = dudx + dvdy;
            Real detadx = gradeta(i,j,0,0);
            Real detady = gradeta(i,j,0,1);
            Real t0 = detadx * div - detadx*dudx - detady*dudy;
            Real t1 = detady * div - detadx*dvdx - detady*dvdy;
            Ax(i,j,0,0) += t0;
            Ax(i,j,0,1) += t1;
        }
    }
}

}

#endif
