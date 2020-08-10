#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
    : amrex::Gpu::Managed
{
    amrex::Real rho_1 = 0.5;
    amrex::Real rho_2 = 2.0;
    amrex::Real p0_base = 5.0;
};

#endif
