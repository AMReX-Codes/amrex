
#include "cns_prob.H"
#include "CNS_K.H"
#include "CNS_index_macros.H"

using namespace amrex;

extern "C" {
    void amrex_probinit (const int* init,
                         const int* name,
                         const int* namelen,
                         const amrex_real* problo,
                         const amrex_real* probhi) {}
}


AMREX_GPU_DEVICE
void
cns_initdata (Box const& bx, FArrayBox& statefab, GeometryData const& geomdata)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto state = statefab.view(lo);
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx      = geomdata.CellSize();

    const Real gamma = 1.4;

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                Real x = prob_lo[0] + (i+lo.x+0.5)*dx[0];
                Real Pt, rhot, uxt;
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
                Real et = Pt/(gamma-1.0);
                state(i,j,k,UEINT) = et;
                state(i,j,k,UEDEN) = et + 0.5*rhot*uxt*uxt;
                state(i,j,k,UTEMP) = 0.0;
            }
        }
    }
}

