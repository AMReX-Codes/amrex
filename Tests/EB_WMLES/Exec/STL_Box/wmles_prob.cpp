#include <AMReX_PROB_AMR_F.H>
#include <AMReX_ParmParse.H>
#include "WMLES_index_macros.H"
#include "WMLES_parm.H"
#include "WMLES.H"

using namespace amrex;

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* /*problo*/,
                         const amrex_real* /*probhi*/)
    {
        ParmParse pp("prob");

        pp.query("inflow_velocity", WMLES::h_prob_parm->inflow_velocity);
        pp.query("inflow_rho",      WMLES::h_prob_parm->inflow_rho);
        pp.query("inflow_T",        WMLES::h_prob_parm->inflow_T);
        pp.query("inflow_p",        WMLES::h_prob_parm->inflow_p);

        Gpu::htod_memcpy_async(WMLES::d_prob_parm, WMLES::h_prob_parm, sizeof(ProbParm));

        Gpu::HostVector<Real> inflow_state(WMLES::numState());

        Real rho  = WMLES::h_prob_parm->inflow_rho;
        Real p    = WMLES::h_prob_parm->inflow_p;
        Real u    = WMLES::h_prob_parm->inflow_velocity;
        Real rhoe = p / (WMLES::h_parm->eos_gamma - 1.0);

        inflow_state[URHO ] = rho;
        inflow_state[UMX  ] = rho * u;
        inflow_state[UMY  ] = 0.0;
        inflow_state[UMZ  ] = 0.0;
        inflow_state[UEDEN] = rhoe + 0.5 * rho * u * u;
        inflow_state[UEINT] = rhoe;
        inflow_state[UTEMP] = WMLES::h_prob_parm->inflow_T;

        Gpu::copyAsync(Gpu::hostToDevice, inflow_state.data(),
                       inflow_state.data() + WMLES::numState(),
                       WMLES::h_prob_parm->inflow_state);
        Gpu::streamSynchronize();
    }
}
