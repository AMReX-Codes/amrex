#include "wmles_prob_parm.H"
#include "WMLES.H"
#include "WMLES_index_macros.H"

#include <AMReX_Arena.H>

ProbParm::ProbParm ()
{
    inflow_state = (amrex::Real*)The_Arena()->alloc(sizeof(Real)*NUM_STATE);
}

ProbParm::~ProbParm ()
{
    The_Arena()->free(inflow_state);
}
