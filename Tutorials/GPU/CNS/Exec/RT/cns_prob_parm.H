#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm
{
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_1;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_2;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real p0_base;
}

#endif
