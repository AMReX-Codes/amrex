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
cns_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state, amrex::GeometryData const& geomdata)
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx      = geomdata.CellSize();

    amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];
    amrex::Real Pt, rhot, uxt;
    if (x < 0.5) {
        Pt = ProbParm::p_l;
        rhot = ProbParm::rho_l;
        uxt = ProbParm::u_l;
    } else {
        Pt = ProbParm::p_r;
        rhot = ProbParm::rho_r;
        uxt = ProbParm::u_r;
    }
    state(i,j,k,URHO ) = rhot;
    state(i,j,k,UMX  ) = rhot*uxt;
    state(i,j,k,UMY  ) = 0.0;
    state(i,j,k,UMZ  ) = 0.0;
    amrex::Real et = Pt/(Parm::eos_gamma-1.0);
    state(i,j,k,UEINT) = et;
    state(i,j,k,UEDEN) = et + 0.5*rhot*uxt*uxt;
    state(i,j,k,UTEMP) = 0.0;
}

#endif
