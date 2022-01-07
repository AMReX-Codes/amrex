#ifndef CNS_PROB_PARM_H_
#define CNS_PROB_PARM_H_

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real inflow_T    = 300.0;
    amrex::Real inflow_p    = 1.0e6;
    amrex::Real inflow_mach = 0.8;
    amrex::Real interior_T  = 1500.0;
    amrex::Real interior_p  = 1.0e6;

    amrex::Real* inflow_state = nullptr;

    ProbParm ();
    ~ProbParm ();
};

#endif
