
#include "cns_prob.H"
#include "CNS_K.H"
#include "CNS_index.H"

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
cns_initdata (Box const& bx, FArrayBox& state, GeometryData const& geomdata)
{
    const auto dp0 = state.stridedPtr(bx);
    const auto len = bx.length3d();
    const auto lo  = bx.loVect3d();
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx      = geomdata.CellSize();

    const Real gamma = 1.4;

    for         (int k = 0; k < len[2]; ++k) {
        for     (int j = 0; j < len[1]; ++j) {
            for (int i = 0; i < len[0]; ++i) {
                Real x = prob_lo[0] + (i+lo[0]+0.5)*dx[0];
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
                Real* AMREX_RESTRICT rho = dp0(i,j,k,URHO);
                Real* AMREX_RESTRICT mx  = dp0(i,j,k,UMX);
                Real* AMREX_RESTRICT my  = dp0(i,j,k,UMY);
                Real* AMREX_RESTRICT mz  = dp0(i,j,k,UMY);
                Real* AMREX_RESTRICT ei  = dp0(i,j,k,UEINT);
                Real* AMREX_RESTRICT ed  = dp0(i,j,k,UEDEN);
                Real* AMREX_RESTRICT T   = dp0(i,j,k,UTEMP);

                *rho = rhot;
                *mx  = rhot*uxt;
                *my  = 0.0;
                *mz  = 0.0;
                *ei  = Pt/(gamma-1.0);
                *ed  = (*ei) + 0.5*rhot*uxt*uxt;
                *T   = 0.0;
            }
        }
    }
}

