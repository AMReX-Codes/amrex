#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include "CNS_index_macros.H"

AMREX_GPU_DEVICE
inline
void
cns_initdata (amrex::Box const& bx, amrex::FArrayBox& statefab, amrex::GeometryData const& geomdata)
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto state = statefab.view(lo);
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx      = geomdata.CellSize();

    const amrex::Real gamma = 1.4;

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                amrex::Real x = prob_lo[0] + (i+lo.x+0.5)*dx[0];
                amrex::Real Pt, rhot, uxt;
                if (x < 0.5) {
                    Pt = 1.0;  // xxxxx parameterize them later
                    rhot = 1.0;
                    uxt = 0.0;
                } else {
                    Pt = 0.1;
                    rhot = 0.125;
                    uxt = 0.0;
                }
                state(i,j,k,URHO ) = rhot;
                state(i,j,k,UMX  ) = rhot*uxt;
                state(i,j,k,UMY  ) = 0.0;
                state(i,j,k,UMZ  ) = 0.0;
                amrex::Real et = Pt/(gamma-1.0);
                state(i,j,k,UEINT) = et;
                state(i,j,k,UEDEN) = et + 0.5*rhot*uxt*uxt;
                state(i,j,k,UTEMP) = 0.0;
            }
        }
    }
}

#endif
