#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real p0 = 1.95e4;
    amrex::Real p1 = 7.42e6;
    amrex::Real rho0 = 3.29e-5;
    amrex::Real rho1 = 3.61e-4;
    amrex::Real v0 = 0.0;
    amrex::Real v1 = 0.0;
    amrex::Real x1 = 4.2;
};

#endif
