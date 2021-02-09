#include <CNS.H>
#include <AMReX_PROB_AMR_F.H>
#include <AMReX_ParmParse.H>

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* /*problo*/,
                         const amrex_real* /*probhi*/)
    {
        // could read parmparse parameters here

        amrex::Gpu::htod_memcpy(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
    }
}
