#ifndef _PROB_H_
#define _PROB_H_

#include <AmrCoreAdv_F.H>

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_INLINE
void
initdata(Box const& bx, Array4<Real> const& phi, GeometryData const& geomdata)
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);    

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize();

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            Real z = prob_lo[2] + (0.5+k) * dx[2];
            Real y = prob_lo[1] + (0.5+j) * dx[1];
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real x = prob_lo[0] + (0.5+i) * dx[0]; 
#if (AMREX_SPACEDIM == 2)
                Real r2 = (pow(x-0.5, 2) + pow((y-0.75),2)) / 0.01;
#else
                Real r2 = (pow(x-0.5, 2) + pow((y-0.75),2) + pow((z-0.5),2)) / 0.01;
#endif
                phi(i,j,k) = 1.0 + std::exp(-r2);
            }
        }
    }
}

#endif
