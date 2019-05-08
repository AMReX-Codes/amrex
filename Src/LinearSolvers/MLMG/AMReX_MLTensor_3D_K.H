#ifndef AMREX_MLTENSOR_3D_K_H_
#define AMREX_MLTENSOR_3D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mltensor_adotx_add_extra (Box const& box, Array4<Real> const& Ax,
                               Array4<Real const> const& vel,
                               Array4<Real const> const& gradeta,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
    const Real dxi = 0.5*dxinv[0];
    const Real dyi = 0.5*dxinv[1];
    const Real dzi = 0.5*dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudx = dxi*(vel(i+1,j  ,k  ,0) - vel(i-1,j  ,k  ,0));
                Real dudy = dyi*(vel(i  ,j+1,k  ,0) - vel(i  ,j-1,k  ,0));
                Real dudz = dzi*(vel(i  ,j  ,k+1,0) - vel(i  ,j  ,k-1,0));
                Real dvdx = dxi*(vel(i+1,j  ,k  ,1) - vel(i-1,j  ,k  ,1));
                Real dvdy = dyi*(vel(i  ,j+1,k  ,1) - vel(i  ,j-1,k  ,1));
                Real dvdz = dzi*(vel(i  ,j  ,k+1,1) - vel(i  ,j  ,k-1,1));
                Real dwdx = dxi*(vel(i+1,j  ,k  ,2) - vel(i-1,j  ,k  ,2));
                Real dwdy = dyi*(vel(i  ,j+1,k  ,2) - vel(i  ,j-1,k  ,2));
                Real dwdz = dzi*(vel(i  ,j  ,k+1,2) - vel(i  ,j  ,k-1,2));
                Real div = dudx + dvdy + dwdz;
                Real detadx = gradeta(i,j,k,0);
                Real detady = gradeta(i,j,k,1);
                Real detadz = gradeta(i,j,k,2);
                Real t0 = detadx * div - detadx*dudx - detady*dudy - detadz*dudz;
                Real t1 = detady * div - detadx*dvdx - detady*dvdy - detadz*dvdz;
                Real t2 = detadz * div - detadx*dwdx - detady*dwdy - detadz*dwdz;
                Ax(i,j,k,0) += t0;
                Ax(i,j,k,1) += t1;
                Ax(i,j,k,2) += t2;
            }
        }
    }
}

}

#endif
