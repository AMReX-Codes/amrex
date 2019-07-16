#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <AMReX_REAL.H>

namespace Params
{
    // This is designed to represent MFiX-like conditions where the grid spacing is 
    //     roughly 2.5 times the particle diameter.  In main.cpp we set grid spacing to 1
    //     so here we set cutoff to diameter = 1/2.5 --> cutoff = 0.2
    static constexpr amrex::Real cutoff = 0.2  ;
    static constexpr amrex::Real min_r  = 1.e-4;
}

#endif
