#ifndef AMREX_MLPOISSON_3D_K_H_
#define AMREX_MLPOISSON_3D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_adotx (Box const& box, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Real dhx, Real dhy, Real dhz)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                y(i,j,k) = dhx * (x(i-1,j,k) - 2.0*x(i,j,k) + x(i+1,j,k))
                     +     dhy * (x(i,j-1,k) - 2.0*x(i,j,k) + x(i,j+1,k))
                     +     dhz * (x(i,j,k-1) - 2.0*x(i,j,k) + x(i,j,k+1));
            }
        }
    }
}

}

#endif
