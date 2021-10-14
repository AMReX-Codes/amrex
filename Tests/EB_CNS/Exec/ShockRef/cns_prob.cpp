
#include <AMReX_PROB_AMR_F.H>
#include <AMReX_ParmParse.H>
#include "cns_prob_parm.H"
#include "CNS.H"

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* /*problo*/,
                         const amrex_real* /*probhi*/)
    {
        amrex::ParmParse pp("prob");

        pp.query("p0", CNS::h_prob_parm->p0);
        pp.query("p1", CNS::h_prob_parm->p1);
        pp.query("rho0", CNS::h_prob_parm->rho0);
        pp.query("rho1", CNS::h_prob_parm->rho1);
        pp.query("v0", CNS::h_prob_parm->v0);
        pp.query("v1", CNS::h_prob_parm->v1);
        pp.query("x1", CNS::h_prob_parm->x1);

        amrex::Gpu::copy(amrex::Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm+1,
                         CNS::d_prob_parm);
    }
}
