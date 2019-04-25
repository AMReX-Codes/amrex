#ifndef AMREX_YAFLUXREGISTER_1D_K_H_
#define AMREX_YAFLUXREGISTER_1D_K_H_

#include <AMReX_FArrayBox.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void yafluxreg_crseadd (Box const& bx, Array4<Real> const& d, Array4<int const> const& flag,
                        FArrayBox const& fxfab, Real dtdx, int nc) noexcept
{
    auto const fx = fxfab.array();

    auto const lo = amrex::lbound(bx);
    auto const hi = amrex::ubound(bx);

    for (int i = lo.x; i <= hi.x; ++i) {
        if (flag(i,0,0) == amrex_yafluxreg_crse_fine_boundary_cell)
        {
            if (flag(i-1,0,0) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,0,0,n) -= dtdx*fx(i,0,0,n);
                }
            } else if (flag(i+1,0,0) == amrex_yafluxreg_fine_cell) {
                for (int n = 0; n < nc; ++n) {
                    d(i,0,0,n) += dtdx*fx(i+1,0,0,n);
                }
            }
        }
    }
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
            const int i = lo.x;
            const int ii = (i+1)*rr.x;
            Real* AMREX_RESTRICT dp = &(d(i,0,0,n));
            Real tmp = -dtdx*f(ii,0,0,n);
#ifdef _OPENMP
#pragma omp atomic
#endif
            *dp += tmp;
        }
        break;
    }
    default:
    {
        for (int n = 0; n < nc; ++n) {
            const int i = lo.x;
            const int ii = i*rr.x;
            Real* AMREX_RESTRICT dp = &(d(i,0,0,n));
            Real tmp = dtdx*f(ii,0,0,n);
#ifdef _OPENMP
#pragma omp atomic
#endif
            *dp += tmp;
        }
    }
    }
}

}
#endif
