#ifndef AMREX_MLMG_3D_K_H_
#define AMREX_MLMG_3D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mlmg_lin_cc_interp_r2 (Box const& bx, Array4<Real> const& ff,
                            Array4<Real const> const& cc, int nc) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < nc; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kc = k/2;
            const int koff = 2*(k-kc*2)-1;
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jc = j/2;
                const int joff = 2*(j-jc*2)-1;
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ic = i/2;
                    const int ioff = 2*(i-ic*2)-1;
                    ff(i,j,k,n) = 0.421875*cc(ic     ,jc     ,kc     ,n)
                        +         0.140625*cc(ic+ioff,jc     ,kc     ,n)
                        +         0.140625*cc(ic     ,jc+joff,kc     ,n)
                        +         0.140625*cc(ic     ,jc     ,kc+koff,n)
                        +         0.046875*cc(ic     ,jc+joff,kc+koff,n)
                        +         0.046875*cc(ic+ioff,jc     ,kc+koff,n)
                        +         0.046875*cc(ic+ioff,jc+joff,kc     ,n)
                        +         0.015625*cc(ic+ioff,jc+joff,kc+koff,n);
                }
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
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kc = k/4;
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jc = j/4;
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ic = i/4;
                    ff(i,j,k,n) = cc(ic,jc,kc,n);
                }
            }
        }
    }
}

}

#endif
