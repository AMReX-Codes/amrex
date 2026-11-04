#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <AMReX_REAL.H>

namespace Params
{
    // This is designed to represent MFiX-like conditions where the grid spacing is
    //     roughly 2.5 times the particle diameter.  In main.cpp we set grid spacing to 1
    //     so here we set cutoff to diameter = 1/2.5 --> cutoff = 0.2
    //     These are functions, not variables, so that they can be used in device
    //     code also if ParticleReal is a class type (AMReX_PRECISION=FP32MP2).
    AMREX_GPU_HOST_DEVICE constexpr amrex::ParticleReal cutoff () { return amrex::ParticleReal(0.2); }
    AMREX_GPU_HOST_DEVICE constexpr amrex::ParticleReal min_r () { return amrex::ParticleReal(1.e-4); }
}

#endif
