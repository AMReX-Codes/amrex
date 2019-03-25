#ifndef AMREX_MLMG_2D_K_H_
#define AMREX_MLMG_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlmg_lin_cc_interp_r2 (Box const& bx, Array4<Real> const& ff,
                            Array4<Real const> const& cc, int nc) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < nc; ++n) {
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jc = j/2;
            const int joff = 2*(j-jc*2)-1;
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ic = i/2;
                const int ioff = 2*(i-ic*2)-1;
                ff(i,j,0,n) = 0.5625*cc(ic     ,jc     ,0,n)
                    +         0.1875*cc(ic+ioff,jc     ,0,n)
                    +         0.1875*cc(ic     ,jc+joff,0,n)
                    +         0.0625*cc(ic+ioff,jc+joff,0,n);

            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlmg_lin_cc_interp_r4 (Box const& bx, Array4<Real> const& ff,
                            Array4<Real const> const& cc, int nc) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < nc; ++n) {
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jc = j/4;
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ic = i/4;
                ff(i,j,0,n) = cc(ic,jc,0,n);
            }
        }
    }
}

}

#endif
