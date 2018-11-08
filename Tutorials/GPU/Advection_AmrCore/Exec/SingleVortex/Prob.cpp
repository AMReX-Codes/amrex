#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

using namespace amrex;

void
initdata(const int level, const int time,
         Box const& bx, FArrayBox& phifab, GeometryData const& geomdata)
{

    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto phi = phifab.view(lo);
    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize();

    const Real gamma = 1.4;

# pragma omp parallel for collapse(2)
    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            Real z = prob_lo[2] + (k+0.5) * dx[2];
            Real y = prob_lo[1] + (j+0.5) * dx[1];
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                Real x = prob_lo[0] + (i+0.5) * dx[0]; 
#if (AMREX_SPACEDIM == 2)
                Real r2 = (pow(x-0.5, 2) + pow((y-0.75),2)) / 0.01;
#else
                Real r2 = (pow(x-0.5, 2) + pow((y-0.75),2) + pow((z-0.5),2)) / 0.01;
#endif
                phi(i,j,k) = 1 + exp(-r2);
            }
        }
    }
}
