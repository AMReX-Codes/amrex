
#include <CNS_parm.H>

void Parm::Initialize ()
{
    constexpr amrex::Real Ru = 8.31451e7_rt;
    cv = Ru / (eos_mu * (eos_gamma-1.0_rt));
    cp = eos_gamma * Ru / (eos_mu * (eos_gamma-1.0_rt));
}
