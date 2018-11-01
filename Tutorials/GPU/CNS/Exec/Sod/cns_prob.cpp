
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

#define IDX(i,n) i+dp0.nstride*n

    for         (int k = 0; k < len[2]; ++k) {
        for     (int j = 0; j < len[1]; ++j) {
            Real* AMREX_RESTRICT u = dp0(0,j,k);
            for (int i = 0; i < len[0]; ++i) {
                Real x = prob_lo[0] + (i+lo[0]+0.5)*dx[0];
                Real Pt, rhot, uxt;
                Real gamma = 1.4;
                if (x < 0.5) {
                    Pt = 1.0;  // xxxxx parameterize them later
                    rhot = 1.0;
                    uxt = 0.0;
                } else {
                    Pt = 0.1;
                    rhot = 0.125;
                    uxt = 0.0;
                }
                u[IDX(i,URHO) ] = rhot;
                u[IDX(i,UMX)  ] = rhot*uxt;
                u[IDX(i,UMY)  ] = 0.0;
                u[IDX(i,UMZ)  ] = 0.0;
                u[IDX(i,UEINT)] = Pt/(gamma-1.0);
                u[IDX(i,UEDEN)] = u[IDX(i,UEINT)] + 0.5*rhot*uxt*uxt;
                u[IDX(i,UTEMP)] = 0.0;
            }
        }
    }

#undef IDX
}

