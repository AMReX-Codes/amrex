#ifndef PROB_PARM_H_
#define PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

using namespace amrex::literals;

// Define the struct to contain problem-specific variables here.
struct ProbParm
{
    amrex::Real adv_vel[AMREX_SPACEDIM] = {AMREX_D_DECL(1.0, 1.0, 1.0)}; // Default values
};

#endif
