#ifndef AMREX_MLABECLAP_1D_K_H_
#define AMREX_MLABECLAP_1D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlabeclap_adotx (Box const& box, Array4<Real> const& y, Array4<Real const> const& x,
                      Array4<Real const> const& a, Array4<Real const> const& bX,
                      GpuArray<Real,AMREX_SPACEDIM> const& dxinv, Real alpha, Real beta)
{
    const Real dhx = beta*dxinv[0]*dxinv[0];

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
       y(i,0,0) = alpha*a(i,0,0)*x(i,0,0)
           - dhx * (bX(i+1,0,0)*(x(i+1,0,0) - x(i  ,0,0))
                  - bX(i  ,0,0)*(x(i  ,0,0) - x(i-1,0,0)));
    }
}

}
#endif
