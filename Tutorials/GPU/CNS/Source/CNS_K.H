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
cns_estdt (amrex::Box const& bx, amrex::FArrayBox const& statefab,
           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto state = statefab.view(lo);
    amrex::Real dt = std::numeric_limits<amrex::Real>::max();

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                amrex::Real rho = state(i,j,k,URHO);
                amrex::Real mx  = state(i,j,k,UMX);
                amrex::Real my  = state(i,j,k,UMY);
                amrex::Real mz  = state(i,j,k,UMY);
                amrex::Real ei  = state(i,j,k,UEINT);
                amrex::Real rhoinv = 1.0/amrex::max(rho,Parm::smallr);
                amrex::Real vx = mx*rhoinv;
                amrex::Real vy = my*rhoinv;
                amrex::Real vz = mz*rhoinv;
                amrex::Real p = amrex::max((Parm::eos_gamma-1.0)*ei, Parm::smallp);
                amrex::Real cs = std::sqrt(Parm::eos_gamma*p*rhoinv);
                amrex::Real dtx = dx[0]/(std::abs(vx)+cs);
                amrex::Real dty = dx[1]/(std::abs(vy)+cs);
                amrex::Real dtz = dx[2]/(std::abs(vz)+cs);
                dt = amrex::min(dt,amrex::min(dtx,amrex::min(dty,dtz)));
            }
        }
    }
    
    return dt;
}


AMREX_GPU_DEVICE
inline
void
cns_compute_temperature (amrex::Box const& bx, amrex::FArrayBox& statefab)
{
    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto u = statefab.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                amrex::Real rhoinv = 1.0/u(i,j,k,URHO);
                amrex::Real mx = u(i,j,k,UMX);
                amrex::Real my = u(i,j,k,UMY);
                amrex::Real mz = u(i,j,k,UMZ);
                u(i,j,k,UEINT) = u(i,j,k,UEDEN) - 0.5 * rhoinv * (mx*mx+my*my+mz*mz);
                u(i,j,k,UTEMP) = rhoinv * u(i,j,k,UEINT) * (1.0/Parm::cv);
            }
        }
    }
}

#endif
