#include <AMReX_PROB_AMR_F.H>
#include <AMReX_ParmParse.H>
#include "CNS_index_macros.H"
#include "CNS_parm.H"
#include "CNS.H"

using namespace amrex;

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* /*problo*/,
                         const amrex_real* /*probhi*/)
    {
        ParmParse pp("prob");

        pp.query("inflow_T"   , CNS::h_prob_parm->inflow_T);
        pp.query("inflow_p"   , CNS::h_prob_parm->inflow_p);
        pp.query("inflow_mach", CNS::h_prob_parm->inflow_mach);
        pp.query("interior_T" , CNS::h_prob_parm->interior_T);
        pp.query("interior_P" , CNS::h_prob_parm->interior_p);

#ifdef AMREX_USE_GPU
        // Cannot use Gpu::copy because ProbParm is not trivailly copyable.
        Gpu::htod_memcpy_async(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
#else
        std::memcpy(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
#endif

        Gpu::HostVector<Real> inflow_state(CNS::numState());

        Real rhoe = CNS::h_prob_parm->inflow_p / (CNS::h_parm->eos_gamma - 1.0);
        Real rho = rhoe/(CNS::h_parm->cv * CNS::h_prob_parm->inflow_T);
        Real cs = std::sqrt(CNS::h_parm->eos_gamma * CNS::h_prob_parm->inflow_p / rho);
        Real v = CNS::h_prob_parm->inflow_mach * cs;
        inflow_state[URHO ] = rho;
        inflow_state[UMX  ] = 0.0;
        inflow_state[UMY  ] = 0.0;
        inflow_state[UMZ  ] = rho*v;
        inflow_state[UEDEN] = rhoe + 0.5*rho*v*v;
        inflow_state[UEINT] = rhoe;
        inflow_state[UTEMP] = CNS::h_prob_parm->inflow_T;

        Gpu::copyAsync(Gpu::hostToDevice, inflow_state.data(),
                       inflow_state.data() + CNS::numState(),
                       CNS::h_prob_parm->inflow_state);
        Gpu::streamSynchronize();
    }
}
