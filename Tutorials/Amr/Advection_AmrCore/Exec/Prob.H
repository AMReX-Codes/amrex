#ifndef PROB_H_
#define PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
initdata(amrex::Box const& bx, amrex::Array4<amrex::Real> const& phi,
         amrex::GeometryData const& geomdata)
{
    using namespace amrex;

    const auto lo = lbound(bx);
    const auto hi = ubound(bx);    

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize();

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (Gpu::notInLaunchRegion())
#endif
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            Real y = prob_lo[1] + (0.5+j) * dx[1];
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real x = prob_lo[0] + (0.5+i) * dx[0]; 
                Real r2 = (std::pow(x-0.5, 2) + std::pow((y-0.75),2)) / 0.01;
                phi(i,j,k) = 1.0 + std::exp(-r2);
            }
        }
    }
}

#endif
