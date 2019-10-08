#ifndef AMREX_MLMG_3D_K_H_
#define AMREX_MLMG_3D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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

#ifdef AMREX_USE_EB
template <int R>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlmg_eb_cc_interp_r (Box const& bx, Array4<Real> const& ff, Array4<Real const> const& cc,
                          Array4<EBCellFlag const> const& flag, int nc) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < nc; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kc = k/R;
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jc = j/R;
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ic = i/R;
                    if (flag(i,j,k).isCovered()) {
                        ff(i,j,k,n) = 0.0;
                    } else {
                        ff(i,j,k,n) = cc(ic,jc,kc,n);
                    }
                }
            }
        }
    }
}
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlmg_lin_nd_interp (int i, int j, int k, int n, Array4<Real> const& fine,
                         Array4<Real const> const& crse) noexcept
{
    int ic = amrex::coarsen(i,2);
    int jc = amrex::coarsen(j,2);
    int kc = amrex::coarsen(k,2);
    bool i_is_odd = (ic*2 != i);
    bool j_is_odd = (jc*2 != j);
    bool k_is_odd = (kc*2 != k);
    if (i_is_odd and j_is_odd and k_is_odd) {
        // Fine node at center of cell
        fine(i,j,k,n) = 0.125*(crse(ic,  jc,  kc,n) + crse(ic,  jc,  kc+1,n) +
                               crse(ic,  jc+1,kc,n) + crse(ic,  jc+1,kc+1,n) +
                               crse(ic+1,jc,  kc,n) + crse(ic+1,jc,  kc+1,n) +
                               crse(ic+1,jc+1,kc,n) + crse(ic+1,jc+1,kc+1,n));
    } else if (j_is_odd and k_is_odd) {
        // Node on a Y-Z face
        fine(i,j,k,n) = 0.25*(crse(ic,  jc,  kc,n) + crse(ic,  jc,  kc+1,n) +
                              crse(ic,  jc+1,kc,n) + crse(ic,  jc+1,kc+1,n));
    } else if (i_is_odd and k_is_odd) {
        // Node on a Z-X face
        fine(i,j,k,n) = 0.25*(crse(ic,  jc,kc,n) + crse(ic,  jc,kc+1,n) +
                              crse(ic+1,jc,kc,n) + crse(ic+1,jc,kc+1,n));
    } else if (i_is_odd and j_is_odd) {
        // Node on a X-Y face
        fine(i,j,k,n) = 0.25*(crse(ic  ,jc,kc,n) + crse(ic  ,jc+1,kc,n) +
                              crse(ic+1,jc,kc,n) + crse(ic+1,jc+1,kc,n));
    } else if (i_is_odd) {
        // Node on X line
        fine(i,j,k,n) = 0.5*(crse(ic,jc,kc,n) + crse(ic+1,jc,kc,n));
    } else if (j_is_odd) {
        // Node on Y line
        fine(i,j,k,n) = 0.5*(crse(ic,jc,kc,n) + crse(ic,jc+1,kc,n));
    } else if (k_is_odd) {
        // Node on Z line
        fine(i,j,k,n) = 0.5*(crse(ic,jc,kc,n) + crse(ic,jc,kc+1,n));
    } else {
        // Node coincident with coarse node
        fine(i,j,k,n) = crse(ic,jc,kc,n);
    }
}

}
#endif
