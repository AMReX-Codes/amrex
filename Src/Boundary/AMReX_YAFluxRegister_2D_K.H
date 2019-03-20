#ifndef AMREX_YAFLUXREGISTER_2D_K_H_
#define AMREX_YAFLUXREGISTER_2D_K_H_

#include <AMReX_FArrayBox.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void yafluxreg_crseadd (Box const& bx, Array4<Real> const& d, Array4<int const> const& flag,
                        FArrayBox const& fxfab, FArrayBox const& fyfab,
                        Real dtdx, Real dtdy, int nc) noexcept
{
    auto const fx = fxfab.array();
    auto const fy = fyfab.array();

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
            } else if (flag(i+1,j,0) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,j,0,n) += dtdx*fx(i+1,j,0,n);
                }
            }

            if (flag(i,j-1,0) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,j,0,n) -= dtdy*fy(i,j,0,n);
                }
            } else if (flag(i,j+1,0) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,j,0,n) += dtdy*fy(i,j+1,0,n);
                }
            }
        }
    }}
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void yafluxreg_fineadd (Box const& bx, Array4<Real> const& d, FArrayBox const& ffab,
                        Real dtdx, int nc, int dirside, Dim3 const& rr) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    const auto f = ffab.array();

    switch (dirside) {
    case 0 :
    {
        for (int n = 0; n < nc; ++n) {
            for (int j = lo.y; j <= hi.y; ++j) {
                const int i = lo.x;
                const int ii = (i+1)*rr.x;
                Real* AMREX_RESTRICT dp = &(d(i,j,0,n));
                for (int joff = 0; joff < rr.y; ++joff) {
                    const int jj =      j*rr.y  + joff;
                    Real tmp = -dtdx*f(ii,jj,0,n);
#ifdef _OPENMP
#pragma omp atomic
#endif
                    *dp += tmp;
                }
            }
        }
        break;
    }
    case 1 :
    {
        for (int n = 0; n < nc; ++n) {
            for (int j = lo.y; j <= hi.y; ++j) {
                const int i = lo.x;
                const int ii = i*rr.x;
                Real* AMREX_RESTRICT dp = &(d(i,j,0,n));
                for (int joff = 0; joff < rr.y; ++joff) {
                    const int jj =      j*rr.y  + joff;
                    Real tmp = dtdx*f(ii,jj,0,n);
#ifdef _OPENMP
#pragma omp atomic
#endif
                    *dp += tmp;
                }
            }
        }
        break;
    }
    case 2 :
    {
        for (int n = 0; n < nc; ++n) {
            for (int i = lo.x; i <= hi.x; ++i) {
                const int j = lo.y;
                const int jj = (j+1)*rr.y;
                Real* AMREX_RESTRICT dp = &(d(i,j,0,n));
                for (int ioff = 0; ioff < rr.x; ++ioff) {
                    const int ii =      i*rr.x  + ioff;
                    Real tmp = -dtdx*f(ii,jj,0,n);
#ifdef _OPENMP
#pragma omp atomic
#endif
                    *dp += tmp;
                }
            }
        }
        break;
    }
    default:
    {
        for (int n = 0; n < nc; ++n) {
            for (int i = lo.x; i <= hi.x; ++i) {
                const int j = lo.y;
                const int jj = j*rr.y;
                Real* AMREX_RESTRICT dp = &(d(i,j,0,n));
                for (int ioff = 0; ioff < rr.x; ++ioff) {
                    const int ii =      i*rr.x  + ioff;
                    Real tmp = dtdx*f(ii,jj,0,n);
#ifdef _OPENMP
#pragma omp atomic
#endif
                    *dp += tmp;
                }
            }
        }
    }
    }
}

}
#endif
