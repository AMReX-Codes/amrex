#ifndef CNS_PARM_H_
#define CNS_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct Parm
    : public amrex::Gpu::Managed
{
    amrex::Real eos_gamma = 1.4_rt;
    amrex::Real eos_mu = 28.97_rt;  // mean molecular weight

    amrex::Real cv;
    amrex::Real cp;

    amrex::Real Pr  = 0.72_rt;     // Prandtl number
    amrex::Real C_S = 1.458e-5_rt; // constant in Sutherland's law
    amrex::Real T_S = 110.4_rt;    // Sutherland temperature

    amrex::Real smallr = 1.e-19_rt;
    amrex::Real smallp = 1.e-10_rt;

    void Initialize ();
};

#endif
