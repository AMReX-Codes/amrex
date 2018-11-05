
#include <limits>
#include <cmath>

#include "CNS_K.H"
#include "CNS_index_macros.H"

using namespace amrex;

AMREX_GPU_DEVICE
Real
cns_estdt (amrex::Box const& bx, amrex::FArrayBox const& statefab,
           GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto state = statefab.view(lo);
    Real dt = std::numeric_limits<Real>::max();
    const Real smallr = 1.e-19;
    const Real smallp = 1.e-10;
    const Real gamma = 1.4;

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                Real rho = state(i,j,k,URHO);
                Real mx  = state(i,j,k,UMX);
                Real my  = state(i,j,k,UMY);
                Real mz  = state(i,j,k,UMY);
                Real ei  = state(i,j,k,UEINT);
                Real rhoinv = 1.0/((rho > smallr) ? rho : smallr);
                Real vx = mx*rhoinv;
                Real vy = my*rhoinv;
                Real vz = mz*rhoinv;
                Real p = (gamma-1.0)*ei;
                p = (p > smallp) ? p : smallp;
                Real cs = std::sqrt(gamma*p*rhoinv);
                Real dtx = dx[0]/(std::abs(vx)+cs);
                Real dty = dx[1]/(std::abs(vy)+cs);
                Real dtz = dx[2]/(std::abs(vz)+cs);
                dt = (dt < dtx) ? dt : dtx;
                dt = (dt < dty) ? dt : dty;
                dt = (dt < dtz) ? dt : dtz;
            }
        }
    }
    
    return dt;
}

