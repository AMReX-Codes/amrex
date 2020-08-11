#ifndef AMREX_MULTIFAB_UTIL_ND_C_H_
#define AMREX_MULTIFAB_UTIL_ND_C_H_

#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <cmath>

namespace amrex {

AMREX_GPU_HOST_DEVICE
inline
void amrex_fill_slice_interp (Box const& bx, Array4<Real> slice,
                              Array4<Real const> const& full,
                              int scomp, int fcomp, int ncomp,
                              int dir, Real coord, GeometryData const& gd) noexcept
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    int ilo = 0, jlo = 0, klo = 0;
    int ihi = 0, jhi = 0, khi = 0;

    Real weight = (coord - gd.ProbLo(dir)) / gd.CellSize(dir);
    int dirhi = static_cast<int>(amrex::Math::floor(weight + Real(0.5)));
    int dirlo = dirhi-1;
    weight -= dirlo+Real(0.5);

    if (weight < 0.5) {
        switch (dir)
        {
        case 0:
            ihi = 1;
            break;
        case 1:
            jhi = 1;
            break;
        case 2:
            khi = 1;
        }
    } else {
        switch (dir)
        {
        case 0:
            ilo = -1;
            break;
        case 1:
            jlo = -1;
            break;
        case 2:
            klo = -1;
        }
    }

    for (int n = 0; n < ncomp; ++n) {
        const int ns = n + scomp;
        const int nf = n + fcomp;
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slice(i,j,k,ns) = (1.0-weight)*full(i+ilo,j+jlo,k+klo,nf)
                        +                  weight *full(i+ihi,j+jhi,k+khi,nf);
                }
            }
        }
    }
}

}

#endif
