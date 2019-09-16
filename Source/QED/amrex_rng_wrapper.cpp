//This file provides a wrapper aroud the RNG
//provided by the amrex library

#include "amrex_rng_wrapper.h"

//RNG wrapper BW engine
//Get rnd number uniformly distributed in [a,b)
amrex::Real AMREX_GPU_DEVICE
amrex_rng_wrapper::unf(amrex::Real a, amrex::Real b)
{
        return (b-a)*amrex::Random() + a;
}

//Get rnd number with exponential distribution
amrex::Real AMREX_GPU_DEVICE
amrex_rng_wrapper::exp(amrex::Real l)
{
    amrex::Real zero_plus_to_one = 1.0 - unf(0.0, 1.0);
    return -log(zero_plus_to_one)/l;
}
