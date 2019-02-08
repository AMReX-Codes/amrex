#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <AMReX_REAL.H>

namespace Params
{
    static constexpr amrex::Real cutoff = 1.e-2;
    static constexpr amrex::Real min_r  = 1.e-4;
    static constexpr amrex::Real mass   = 1.e-2;
}

#endif
