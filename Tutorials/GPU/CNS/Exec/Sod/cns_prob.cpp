
#include <AMReX_PROB_AMR_F.H>
#include <AMReX_ParmParse.H>
#include "cns_prob_parm.H"

namespace ProbParm
{
    AMREX_GPU_DEVICE_MANAGED amrex::Real p_l = 1.0;
    AMREX_GPU_DEVICE_MANAGED amrex::Real p_r = 0.1;
    AMREX_GPU_DEVICE_MANAGED amrex::Real rho_l = 1.0;
    AMREX_GPU_DEVICE_MANAGED amrex::Real rho_r = 0.125;
    AMREX_GPU_DEVICE_MANAGED amrex::Real u_l = 0.0;
    AMREX_GPU_DEVICE_MANAGED amrex::Real u_r = 0.0;
}

extern "C" {
    void amrex_probinit (const int* init,
                         const int* name,
                         const int* namelen,
                         const amrex_real* problo,
                         const amrex_real* probhi)
    {
        amrex::ParmParse pp("prob");

        pp.query("p_l", ProbParm::p_l);
        pp.query("p_r", ProbParm::p_r);
        pp.query("rho_l", ProbParm::rho_l);
        pp.query("rho_r", ProbParm::rho_r);
        pp.query("u_l", ProbParm::u_l);
        pp.query("u_r", ProbParm::u_r);
    }
}
