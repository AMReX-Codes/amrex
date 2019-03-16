#ifndef AMREX_MLMG_1D_K_H_
#define AMREX_MLMG_1D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlmg_lin_cc_interp_r2 (Box const& bx, Array4<Real> const& ff,
                            Array4<Real const> const& cc, int nc) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < nc; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            const int ic = i/2;
            const int ioff = 2*(i-ic*2)-1;
            ff(i,0,0,n) = 0.75*cc(ic,0,0,n) + 0.25*cc(ic+ioff,0,0,n);
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
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            const int ic = i/4;
            ff(i,0,0,n) = cc(ic,0,0,n);
        }
    }
}

}

#endif
