
#include <CNS.H>

using namespace amrex;

Real
CNS::advance (Real time, Real dt, int iteration, int ncycle)
{
    BL_PROFILE("CNS::advance()");
}

void
CNS::compute_dSdt (const MultiFab& S, MultiFab& dSdt, Real dt,
                   YAFluxRegister* fr_as_crse, YAFluxRegister* fr_as_fine)
{
    BL_PROFILE("CNS::compute_dSdt()");
}
