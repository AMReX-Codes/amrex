#ifndef AMREX_MLMG_2D_K_H_
#define AMREX_MLMG_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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

#ifdef AMREX_USE_EB
template <int R>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlmg_eb_cc_interp_r (Box const& bx, Array4<Real> const& ff, Array4<Real const> const& cc,
                          Array4<EBCellFlag const> const& flag, int nc) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < nc; ++n) {
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jc = j/R;
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ic = i/R;
                if (flag(i,j,0).isCovered()) {
                    ff(i,j,0,n) = 0.0;
                } else {
                    ff(i,j,0,n) = cc(ic,jc,0,n);
                }
            }
        }
    }
}
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlmg_lin_nd_interp (int i, int j, int, int n, Array4<Real> const& fine,
                         Array4<Real const> const& crse) noexcept
{
    int ic = amrex::coarsen(i,2);
    int jc = amrex::coarsen(j,2);
    bool i_is_odd = (ic*2 != i);
    bool j_is_odd = (jc*2 != j);
    if (i_is_odd and j_is_odd) {
        // Fine node at center of cell
        fine(i,j,0,n) = 0.25*(crse(ic,jc,0,n) + crse(ic+1,jc,0,n) +
                              crse(ic,jc+1,0,n) + crse(ic+1,jc+1,0,n));
    } else if (j_is_odd) {
        // Fine node along a vertical cell edge
        fine(i,j,0,n) = 0.5*(crse(ic,jc,0,n) + crse(ic,jc+1,0,n));
    } else if (i_is_odd) {
        // Fine node along a horizontal cell edge
        fine(i,j,0,n) = 0.5*(crse(ic,jc,0,n) + crse(ic+1,jc,0,n));
    } else {
        // Fine node coincident with coarse node
        fine(i,j,0,n) = crse(ic,jc,0,n);
    }
}

}

#endif
