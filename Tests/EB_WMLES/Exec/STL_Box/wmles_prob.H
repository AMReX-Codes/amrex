#ifndef WMLES_PROB_H_
#define WMLES_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include "WMLES_index_macros.H"
#include "WMLES_parm.H"
#include "wmles_prob_parm.H"

AMREX_GPU_DEVICE
inline
void
wmles_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& prob_parm)
{
    using amrex::Real;

    Real rho   = prob_parm.inflow_rho;
    Real u_inf = prob_parm.inflow_velocity;

    state(i,j,k,URHO)  = rho;
    state(i,j,k,UMX)   = rho * u_inf;
    state(i,j,k,UMY)   = 0.0;
    state(i,j,k,UMZ)   = 0.0;

    Real p   = prob_parm.inflow_p;
    Real rhoe = p / (parm.eos_gamma - Real(1.0));

    state(i,j,k,UEINT) = rhoe;
    state(i,j,k,UEDEN) = rhoe + Real(0.5) * rho * u_inf * u_inf;
    state(i,j,k,UTEMP) = prob_parm.inflow_T;
}

#endif
