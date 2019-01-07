
#include <AMReX_PROB_AMR_F.H>
#include "cns_prob_parm.H"
#include <AMReX_ParmParse.H>

namespace ProbParm
{
    AMREX_GPU_DEVICE_MANAGED amrex::Real rho_1 = 0.5;
    AMREX_GPU_DEVICE_MANAGED amrex::Real rho_2 = 2.0;
    AMREX_GPU_DEVICE_MANAGED amrex::Real p0_base = 5.0;
}

extern "C" {
    void amrex_probinit (const int* init,
                         const int* name,
                         const int* namelen,
                         const amrex_real* problo,
                         const amrex_real* probhi)
    {
        // could read parmparse parameters here
    }
}
