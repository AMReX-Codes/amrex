#ifndef AMREX_YAFLUXREGISTER_2D_K_H_
#define AMREX_YAFLUXREGISTER_2D_K_H_
#include <AMReX_Config.H>

#include <AMReX_BaseFab.H>

namespace amrex {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void yafluxreg_crseadd (Box const& bx, Array4<T> const& d, Array4<int const> const& flag,
                        Array4<T const> const& fx, Array4<T const> const& fy,
                        T dtdx, T dtdy, int nc) noexcept
{
    auto const lo = amrex::lbound(bx);
    auto const hi = amrex::ubound(bx);

    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        if (flag(i,j,0) == amrex_yafluxreg_crse_fine_boundary_cell)
        {
            if (flag(i-1,j,0) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,j,0,n) -= dtdx*fx(i,j,0,n);
                }
            }
            if (flag(i+1,j,0) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,j,0,n) += dtdx*fx(i+1,j,0,n);
                }
            }

            if (flag(i,j-1,0) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,j,0,n) -= dtdy*fy(i,j,0,n);
                }
            }
            if (flag(i,j+1,0) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,j,0,n) += dtdy*fy(i,j+1,0,n);
                }
            }
        }
    }}
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T yafluxreg_fineadd_delta (int i, int j, int n,
                           Array4<T const> const& f, T scale,
                           int dirside, Dim3 const& rr) noexcept
{
    const int dir  = dirside / 2;
    const int side = dirside % 2;

    T sum = T(0.0);

    if (dir == 0) {
        const int ii    = (side == 0) ? (i+1) * rr.x : i * rr.x;
        const T sign    = (side == 0) ? T(-1.0) : T(1.0);
        const int jbase = j * rr.y;
#pragma unroll
        for (int jj = 0; jj < rr.y; ++jj) {
            const int fj = jbase + jj;
            sum += sign * scale * f(ii, fj, 0, n);
        }
    } else {
        const int jj    = (side == 0) ? (j+1) * rr.y : j * rr.y;
        const T sign    = (side == 0) ? T(-1.0) : T(1.0);
        const int ibase = i * rr.x;
#pragma unroll
        for (int ii = 0; ii < rr.x; ++ii) {
            const int fi = ibase + ii;
            sum += sign * scale * f(fi, jj, 0, n);
        }
    }

    return sum;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void yafluxreg_fineadd (int i, int j, int n,
                        Array4<T> const& d, Array4<T const> const& f,
                        T scale, int dirside, Dim3 const& rr) noexcept
{
    d(i,j,0,n) += yafluxreg_fineadd_delta(i,j,n,f,scale,dirside,rr);
}

}
#endif
