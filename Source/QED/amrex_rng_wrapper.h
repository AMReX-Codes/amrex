#ifndef WARPX_amrex_rng_wrapper_h_
#define WARPX_amrex_rng_wrapper_h_

//This file provides a wrapper aroud the RNG
//provided by the amrex library

#include <AMReX_AmrCore.H>

//RNG wrapper BW engine
class amrex_rng_wrapper
{
public:
    //Get rnd number uniformly distributed in [a,b)
    amrex::Real AMREX_GPU_DEVICE
    unf(amrex::Real a, amrex::Real b);

    //Get rnd number with exponential distribution
    amrex::Real AMREX_GPU_DEVICE
    exp(amrex::Real l);
};

#endif //WARPX_amrex_rng_wrapper_h_
