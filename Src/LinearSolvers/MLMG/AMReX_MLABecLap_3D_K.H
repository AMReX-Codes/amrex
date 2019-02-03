#ifndef AMREX_MLABECLAP_3D_K_H_
#define AMREX_MLABECLAP_3D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_adotx (Box const& box, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Array4<Real const> const& a,
                      Array4<Real const> const& bX,
                      Array4<Real const> const& bY,
                      Array4<Real const> const& bZ,
                      GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                      Real alpha, Real beta)
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];
    const Real dhz = beta*dxinv[2]*dxinv[2];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                y(i,j,k) = alpha*a(i,j,k)*x(i,j,k)
                    - dhx * (bX(i+1,j,k)*(x(i+1,j,k) - x(i  ,j,k))
                           - bX(i  ,j,k)*(x(i  ,j,k) - x(i-1,j,k)))
                    - dhy * (bY(i,j+1,k)*(x(i,j+1,k) - x(i,j  ,k))
                           - bY(i,j  ,k)*(x(i,j  ,k) - x(i,j-1,k)))
                    - dhz * (bZ(i,j,k+1)*(x(i,j,k+1) - x(i,j,k  ))
                           - bZ(i,j,k  )*(x(i,j,k  ) - x(i,j,k-1)));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_normalize (Box const& box, Array4<Real> const& x,
                          Array4<Real const> const& a,
                          Array4<Real const> const& bX,
                          Array4<Real const> const& bY,
                          Array4<Real const> const& bZ,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                          Real alpha, Real beta)
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];
    const Real dhz = beta*dxinv[2]*dxinv[2];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                x(i,j,k) /= alpha*a(i,j,k)
                    + dhx*(bX(i,j,k)+bX(i+1,j,k))
                    + dhy*(bY(i,j,k)+bY(i,j+1,k))
                    + dhz*(bZ(i,j,k)+bZ(i,j,k+1));
            }
        }
    }
}

}
#endif
