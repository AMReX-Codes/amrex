#ifndef WMLES_PROB_PARM_H_
#define WMLES_PROB_PARM_H_

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real inflow_velocity = 1.0;
    amrex::Real inflow_rho     = 1.0;
    amrex::Real inflow_T       = 300.0;
    amrex::Real inflow_p       = 1.0e5;

    amrex::Real* inflow_state = nullptr;

    ProbParm ();
    ~ProbParm ();
};

#endif
