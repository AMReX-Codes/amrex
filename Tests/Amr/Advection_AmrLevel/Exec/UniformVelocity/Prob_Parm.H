#ifndef PROB_PARM_H_
#define PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real adv_vel[3] = {1.0, 1.0, 1.0};
};

ProbParm* h_prob_parm;
ProbParm* d_prob_parm;

#endif
