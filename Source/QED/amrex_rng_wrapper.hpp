#ifndef WARPX_amrex_rng_wrapper_hpp_
#define WARPX_rng_wrapper_hpp_

//This file provides a wrapper aroud the RNG
//provided by the amrex library

#include <AMReX_AmrCore.H>

//RNG wrapper BW engine
inline class amrex_rng_wrapper
{
public:
    //Get rnd number uniformly distributed in [a,b)
    amrex::Real unf(amrex::Real a, amrex::Real b)
    {
        return (b-a)*amrex::Random() + a;
    }

    //Get rnd number with exponential distribution
    amrex::Real exp(amrex::Real l)
    {
        amrex::Real zero_plus_to_one = 1.0 - unf(0.0, 1.0);
        return -log(zero_plus_to_one)/l;
    }
};

#endif //amrex_rng_wrapper_hpp_
