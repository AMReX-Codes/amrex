#ifndef CNS_PARM_H_
#define CNS_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace Parm
{
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real eos_gamma;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real eos_mu;

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real cv;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real cp;

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real Pr;  // Prandtl number
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real C_S; // constant in Sutherland's law
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real T_S; // Sutherland temperature

    extern AMREX_GPU_DEVICE_MANAGED amrex::Real smallr;
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real smallp;

    void Initialize ();
}

#endif
