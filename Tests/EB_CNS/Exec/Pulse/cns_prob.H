#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include "CNS_index_macros.H"
#include "CNS_parm.H"
#include "cns_prob_parm.H"

AMREX_GPU_DEVICE
inline
void
cns_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& prob_parm)
{
    using amrex::Real;

    constexpr Real pi = 3.1415926535897932;

    const Real* prob_lo = geomdata.ProbLo();
    const Real* prob_hi = geomdata.ProbHi();
    const Real* dx      = geomdata.CellSize();

    Real x = (prob_lo[0] + (i+Real(0.5))*dx[0]) - 0.5 * (prob_lo[0] + prob_hi[0]);
    Real y = (prob_lo[1] + (j+Real(0.5))*dx[1]) - 0.5 * (prob_lo[1] + prob_hi[1]);
#if (AMREX_SPACEDIM == 3)
    Real z = (prob_lo[2] + (k+Real(0.5))*dx[2]) - 0.5 * (prob_lo[2] + prob_hi[2]);
#else
    Real z = Real(0.0);
#endif
    Real r = std::sqrt(x*x + y*y + z*z);

    if (r > prob_parm.rpulse) {
        state(i,j,k,URHO ) = prob_parm.rho0;
    } else {
        state(i,j,k,URHO ) = prob_parm.rho0 + prob_parm.drho0 * exp(-16.*r*r) * std::pow(std::cos(pi*r),6.0);
    }

    state(i,j,k,UMX  ) = Real(0.0);
    state(i,j,k,UMY  ) = Real(0.0);
#if (AMREX_SPACEDIM == 3)
    state(i,j,k,UMZ  ) = Real(0.0);
#endif

    Real Pt = prob_parm.p0 * std::pow( (state(i,j,k,URHO) / prob_parm.rho0), parm.eos_gamma );

    state(i,j,k,UEINT) = Pt / (parm.eos_gamma - 1.0);
    state(i,j,k,UEDEN) = state(i,j,k,UEINT);
    state(i,j,k,UTEMP) = state(i,j,k,UEINT) / (parm.cv * state(i,j,k,URHO));
}

#endif
