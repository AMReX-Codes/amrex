#ifndef CNS_PARM_H_
#define CNS_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace Parm
{
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real eos_gamma;

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real smallr;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real smallp;

}

#endif
