#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
    : amrex::Gpu::Managed
{
    amrex::Real p_l = 1.0;
    amrex::Real p_r = 0.1;
    amrex::Real rho_l = 1.0;
    amrex::Real rho_r = 0.125;
    amrex::Real u_l = 0.0;
    amrex::Real u_r = 0.0;
};

#endif
