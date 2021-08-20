#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real rpulse = 0.5;
    amrex::Real rho0   = 1.2e-3;
    amrex::Real drho0  = 1.2e-4;
    amrex::Real p0     = 1.01325e6;
};

#endif
