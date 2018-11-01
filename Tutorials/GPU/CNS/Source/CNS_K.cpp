
#include <limits>
#include <cmath>

#include "CNS_K.H"
#include "CNS_index.H"

using namespace amrex;

AMREX_GPU_DEVICE
Real
cns_estdt (amrex::Box const& bx, amrex::FArrayBox const& state,
           GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    const auto dp0 = state.stridedPtr(bx);
    const auto len = bx.length3d();
    Real dt = std::numeric_limits<Real>::max();
    const Real smallr = 1.e-19;
    const Real smallp = 1.e-10;
    const Real gamma = 1.4;

    for         (int k = 0; k < len[2]; ++k) {
        for     (int j = 0; j < len[1]; ++j) {
            for (int i = 0; i < len[0]; ++i) {
                Real rho = *dp0(i,j,k,URHO);
                Real mx  = *dp0(i,j,k,UMX);
                Real my  = *dp0(i,j,k,UMY);
                Real mz  = *dp0(i,j,k,UMY);
                Real ei  = *dp0(i,j,k,UEINT);
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

