#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm
{
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real p_l;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real p_r;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_l;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_r;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real u_l;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real u_r;
}

#endif
