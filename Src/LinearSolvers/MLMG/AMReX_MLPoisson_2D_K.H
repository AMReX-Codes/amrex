#ifndef AMREX_MLPOISSON_2D_K_H_
#define AMREX_MLPOISSON_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlpoisson_adotx (Box const& box, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Real dhx, Real dhy,
                      Real const* AMREX_RESTRICT rc,
                      Real const* AMREX_RESTRICT re,
                      int rlo)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            y(i,j,0) = dhx * (re[i-rlo]*x(i-1,j,0)
                              - (re[i-rlo]+re[i+1-rlo])*x(i,j,0)
                              + re[i+1-rlo]*x(i+1,j,0))
                +   dhy * rc[i-rlo] * (x(i,j-1,0) - 2.*x(i,j,0) + x(i,j+1,0));
        }
    }
}

}

#endif
