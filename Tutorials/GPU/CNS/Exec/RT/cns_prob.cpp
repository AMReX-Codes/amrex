
#include <AMReX_PROB_AMR_F.H>
#include "cns_prob_parm.H"
#include <AMReX_ParmParse.H>

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* /*problo*/,
                         const amrex_real* /*probhi*/)
    {
        // could read parmparse parameters here
    }
}
