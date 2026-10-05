#ifndef AMREX_MLMG_3D_K_H_
#define AMREX_MLMG_3D_K_H_
#include <AMReX_Config.H>

namespace amrex {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlmg_lin_cc_interp_r2 (Box const& bx, Array4<T> const& ff,
                            Array4<T const> const& cc, int nc) noexcept
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
                    ff(i,j,k,n) = T(0.421875)*cc(ic     ,jc     ,kc     ,n)
                        +         T(0.140625)*cc(ic+ioff,jc     ,kc     ,n)
                        +         T(0.140625)*cc(ic     ,jc+joff,kc     ,n)
                        +         T(0.140625)*cc(ic     ,jc     ,kc+koff,n)
                        +         T(0.046875)*cc(ic     ,jc+joff,kc+koff,n)
                        +         T(0.046875)*cc(ic+ioff,jc     ,kc+koff,n)
                        +         T(0.046875)*cc(ic+ioff,jc+joff,kc     ,n)
                        +         T(0.015625)*cc(ic+ioff,jc+joff,kc+koff,n);
                }
            }
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlmg_lin_cc_interp_r4 (Box const& bx, Array4<T> const& ff,
                            Array4<T const> const& cc, int nc) noexcept
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
template <int R, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlmg_eb_cc_interp_r (Box const& bx, Array4<T> const& ff, Array4<T const> const& cc,
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
                        ff(i,j,k,n) = T(0.0);
                    } else {
                        ff(i,j,k,n) = cc(ic,jc,kc,n);
                    }
                }
            }
        }
    }
}
#endif

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlmg_lin_nd_interp_r2 (int i, int j, int k, int n, Array4<T> const& fine,
                            Array4<T const> const& crse) noexcept
{
    int ic = amrex::coarsen(i,2);
    int jc = amrex::coarsen(j,2);
    int kc = amrex::coarsen(k,2);
    bool i_is_odd = (ic*2 != i);
    bool j_is_odd = (jc*2 != j);
    bool k_is_odd = (kc*2 != k);
    if (i_is_odd && j_is_odd && k_is_odd) {
        // Fine node at center of cell
        fine(i,j,k,n) = T(0.125)*(crse(ic,  jc,  kc,n) + crse(ic,  jc,  kc+1,n) +
                                  crse(ic,  jc+1,kc,n) + crse(ic,  jc+1,kc+1,n) +
                                  crse(ic+1,jc,  kc,n) + crse(ic+1,jc,  kc+1,n) +
                                  crse(ic+1,jc+1,kc,n) + crse(ic+1,jc+1,kc+1,n));
    } else if (j_is_odd && k_is_odd) {
        // Node on a Y-Z face
        fine(i,j,k,n) = T(0.25)*(crse(ic,  jc,  kc,n) + crse(ic,  jc,  kc+1,n) +
                                 crse(ic,  jc+1,kc,n) + crse(ic,  jc+1,kc+1,n));
    } else if (i_is_odd && k_is_odd) {
        // Node on a Z-X face
        fine(i,j,k,n) = T(0.25)*(crse(ic,  jc,kc,n) + crse(ic,  jc,kc+1,n) +
                                 crse(ic+1,jc,kc,n) + crse(ic+1,jc,kc+1,n));
    } else if (i_is_odd && j_is_odd) {
        // Node on a X-Y face
        fine(i,j,k,n) = T(0.25)*(crse(ic  ,jc,kc,n) + crse(ic  ,jc+1,kc,n) +
                                 crse(ic+1,jc,kc,n) + crse(ic+1,jc+1,kc,n));
    } else if (i_is_odd) {
        // Node on X line
        fine(i,j,k,n) = T(0.5)*(crse(ic,jc,kc,n) + crse(ic+1,jc,kc,n));
    } else if (j_is_odd) {
        // Node on Y line
        fine(i,j,k,n) = T(0.5)*(crse(ic,jc,kc,n) + crse(ic,jc+1,kc,n));
    } else if (k_is_odd) {
        // Node on Z line
        fine(i,j,k,n) = T(0.5)*(crse(ic,jc,kc,n) + crse(ic,jc,kc+1,n));
    } else {
        // Node coincident with coarse node
        fine(i,j,k,n) = crse(ic,jc,kc,n);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlmg_lin_nd_interp_r4 (int i, int j, int k, int n, Array4<T> const& fine,
                            Array4<T const> const& crse) noexcept
{
    int ic = amrex::coarsen(i,4);
    int jc = amrex::coarsen(j,4);
    int kc = amrex::coarsen(k,4);
    bool i_injection = (ic*4 == i);
    bool j_injection = (jc*4 == j);
    bool k_injection = (kc*4 == k);

#define I_LO (4*(ic+1)-i)
#define J_LO (4*(jc+1)-j)
#define K_LO (4*(kc+1)-k)
#define I_HI (i-4*ic)
#define J_HI (j-4*jc)
#define K_HI (k-4*kc)

    if (i_injection && j_injection && k_injection)
    {
        fine(i,j,k,n) = crse(ic,jc,kc,n);
    }
    else if (i_injection && j_injection)
    {
        fine(i,j,k,n) = T(0.25)*(crse(ic,jc,kc  ,n)*T(K_LO)
                               + crse(ic,jc,kc+1,n)*T(K_HI));
    }
    else if (i_injection && k_injection)
    {
        fine(i,j,k,n) = T(0.25)*(crse(ic,jc  ,kc,n)*T(J_LO)
                               + crse(ic,jc+1,kc,n)*T(J_HI));
    }
    else if (j_injection && k_injection)
    {
        fine(i,j,k,n) = T(0.25)*(crse(ic  ,jc,kc,n)*T(I_LO)
                               + crse(ic+1,jc,kc,n)*T(I_HI));
    }
    else if (i_injection)
    {
        fine(i,j,k,n) = T(0.0625)*(crse(ic,jc  ,kc  ,n)*T(J_LO*K_LO)
                                 + crse(ic,jc+1,kc  ,n)*T(J_HI*K_LO)
                                 + crse(ic,jc  ,kc+1,n)*T(J_LO*K_HI)
                                 + crse(ic,jc+1,kc+1,n)*T(J_HI*K_HI));
    }
    else if (j_injection)
    {
        fine(i,j,k,n) = T(0.0625)*(crse(ic  ,jc,kc  ,n)*T(I_LO*K_LO)
                                 + crse(ic+1,jc,kc  ,n)*T(I_HI*K_LO)
                                 + crse(ic  ,jc,kc+1,n)*T(I_LO*K_HI)
                                 + crse(ic+1,jc,kc+1,n)*T(I_HI*K_HI));

    } else if (k_injection) {
        fine(i,j,k,n) = T(0.0625)*(crse(ic  ,jc  ,kc,n)*T(I_LO*J_LO)
                                 + crse(ic+1,jc  ,kc,n)*T(I_HI*J_LO)
                                 + crse(ic  ,jc+1,kc,n)*T(I_LO*J_HI)
                                 + crse(ic+1,jc+1,kc,n)*T(I_HI*J_HI));
    }
    else
    {
        fine(i,j,k,n) = T(0.015625)*(crse(ic  ,jc  ,kc  ,n)*T(I_LO*J_LO*K_LO)
                                   + crse(ic+1,jc  ,kc  ,n)*T(I_HI*J_LO*K_LO)
                                   + crse(ic  ,jc+1,kc  ,n)*T(I_LO*J_HI*K_LO)
                                   + crse(ic+1,jc+1,kc  ,n)*T(I_HI*J_HI*K_LO)
                                   + crse(ic  ,jc  ,kc+1,n)*T(I_LO*J_LO*K_HI)
                                   + crse(ic+1,jc  ,kc+1,n)*T(I_HI*J_LO*K_HI)
                                   + crse(ic  ,jc+1,kc+1,n)*T(I_LO*J_HI*K_HI)
                                   + crse(ic+1,jc+1,kc+1,n)*T(I_HI*J_HI*K_HI));
    }

#undef I_LO
#undef J_LO
#undef K_LO
#undef I_HI
#undef J_HI
#undef K_HI
}

}
#endif
