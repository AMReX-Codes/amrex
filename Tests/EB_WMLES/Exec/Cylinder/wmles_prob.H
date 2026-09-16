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

    Real x = geomdata.ProbLo(0) + (i + 0.5) * geomdata.CellSize(0);
    Real y = geomdata.ProbLo(1) + (j + 0.5) * geomdata.CellSize(1);
#if (AMREX_SPACEDIM == 3)
    Real z = geomdata.ProbLo(2) + (k + 0.5) * geomdata.CellSize(2);
#endif

    Real rho   = prob_parm.inflow_rho;
    Real u_inf = prob_parm.inflow_velocity;

    state(i,j,k,URHO)  = rho;
    state(i,j,k,UMX)   = rho * u_inf;
    state(i,j,k,UMY)   = 0.0;
#if (AMREX_SPACEDIM == 3)
    state(i,j,k,UMZ)   = 0.0;
#endif

    Real p   = prob_parm.inflow_p;
    Real rhoe = p / (parm.eos_gamma - Real(1.0));

    state(i,j,k,UEINT) = rhoe;
    state(i,j,k,UEDEN) = rhoe + Real(0.5) * rho * u_inf * u_inf;
    state(i,j,k,UTEMP) = prob_parm.inflow_T;
}

#endif
