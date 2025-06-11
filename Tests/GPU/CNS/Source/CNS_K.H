#ifndef CNS_K_H_
#define CNS_K_H_

#include "CNS_index_macros.H"
#include <AMReX_FArrayBox.H>
#include <limits>
#include <cmath>

#include "cns_prob.H"
#include "CNS_parm.H"

AMREX_GPU_HOST_DEVICE
inline
amrex::Real
cns_estdt (amrex::Box const& bx, amrex::Array4<Real const> const& state,
           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
           Parm const& parm) noexcept
{
    using amrex::Real;

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    Real dt = std::numeric_limits<Real>::max();

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real rho = state(i,j,k,URHO);
                Real mx  = state(i,j,k,UMX);
                Real my  = state(i,j,k,UMY);
                Real mz  = state(i,j,k,UMY);
                Real ei  = state(i,j,k,UEINT);
                Real rhoinv = Real(1.0)/amrex::max(rho,parm.smallr);
                Real vx = mx*rhoinv;
                Real vy = my*rhoinv;
                Real vz = mz*rhoinv;
                Real p = amrex::max((parm.eos_gamma-Real(1.0))*ei, parm.smallp);
                Real cs = std::sqrt(parm.eos_gamma*p*rhoinv);
                Real dtx = dx[0]/(std::abs(vx)+cs);
                Real dty = dx[1]/(std::abs(vy)+cs);
                Real dtz = dx[2]/(std::abs(vz)+cs);
                dt = amrex::min(dt,amrex::min(dtx,amrex::min(dty,dtz)));
            }
        }
    }

    return dt;
}


AMREX_GPU_DEVICE
inline
void
cns_compute_temperature (int i, int j, int k, amrex::Array4<amrex::Real> const& u,
                         Parm const& parm) noexcept
{
    using amrex::Real;

    Real rhoinv = Real(1.0)/u(i,j,k,URHO);
    Real mx = u(i,j,k,UMX);
    Real my = u(i,j,k,UMY);
    Real mz = u(i,j,k,UMZ);
    u(i,j,k,UEINT) = u(i,j,k,UEDEN) - Real(0.5) * rhoinv * (mx*mx+my*my+mz*mz);
    u(i,j,k,UTEMP) = rhoinv * u(i,j,k,UEINT) * (Real(1.0)/parm.cv);
}

#endif
