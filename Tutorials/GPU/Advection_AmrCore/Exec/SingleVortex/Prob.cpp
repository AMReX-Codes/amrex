
#include <AmrCoreAdv_F.H>

#include <Prob.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

using namespace amrex;

AMREX_GPU_DEVICE
void
initdata(Box const& bx, FArrayBox& phifab, GeometryData const& geom)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto phi = phifab.view(lo);
    const Real* AMREX_RESTRICT prob_lo = geom.ProbLo();
    const Real* AMREX_RESTRICT dx      = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if (GPU::notInLaunchRegion)
#endif
    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            Real z = prob_lo[2] + (0.5+(k+lo.z)) * dx[2];
            Real y = prob_lo[1] + (0.5+(j+lo.y)) * dx[1];
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                Real x = prob_lo[0] + (0.5+(i+lo.x)) * dx[0]; 
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

