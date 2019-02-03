#ifndef AMREX_MLABECLAP_2D_K_H_
#define AMREX_MLABECLAP_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_adotx (Box const& box, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Array4<Real const> const& a,
                      Array4<Real const> const& bX,
                      Array4<Real const> const& bY,
                      GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                      Real alpha, Real beta)
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            y(i,j,0) = alpha*a(i,j,0)*x(i,j,0)
                - dhx * (bX(i+1,j,0)*(x(i+1,j,0) - x(i  ,j,0))
                       - bX(i  ,j,0)*(x(i  ,j,0) - x(i-1,j,0)))
                - dhy * (bY(i,j+1,0)*(x(i,j+1,0) - x(i,j  ,0))
                       - bY(i,j  ,0)*(x(i,j  ,0) - x(i,j-1,0)));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_normalize (Box const& box, Array4<Real> const& x,
                          Array4<Real const> const& a,
                          Array4<Real const> const& bX,
                          Array4<Real const> const& bY,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                          Real alpha, Real beta)
{
    const Real dhx = beta*dxinv[0]*dxinv[0];
    const Real dhy = beta*dxinv[1]*dxinv[1];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            x(i,j,0) /= alpha*a(i,j,0)
                + dhx*(bX(i,j,0)+bX(i+1,j,0))
                + dhy*(bY(i,j,0)+bY(i,j+1,0));
        }
    }
}

}
#endif
