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
void amrex_fill_slice_interp (Box const& bx, FArrayBox & slicefab, FArrayBox const& fullfab,
                              int scomp, int fcomp, int ncomp,
                              int dir, Real coord, GeometryData const& gd) noexcept
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto full = fullfab.view(lo,fcomp);
    const auto slice = slicefab.view(lo,scomp);

    int ilo = 0, jlo = 0, klo = 0;
    int ihi = 0, jhi = 0, khi = 0;

    Real weight = (coord - gd.ProbLo(dir)) / gd.CellSize(dir);
    int dirhi = std::floor(weight + 0.5);
    int dirlo = dirhi-1;
    weight -= dirlo+0.5;

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
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    slice(i,j,k,n) = (1.0-weight)*full(i+ilo,j+jlo,k+klo,n)
                        +                 weight *full(i+ihi,j+jhi,k+khi,n);
                }
            }
        }
    }
}

}

#endif
