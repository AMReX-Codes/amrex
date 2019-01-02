
#include <CNS_parm.H>

namespace Parm
{
    AMREX_GPU_DEVICE_MANAGED amrex::Real eos_gamma = 1.4;
    AMREX_GPU_DEVICE_MANAGED amrex::Real eos_mu = 28.97;  // mean molecular weight

    AMREX_GPU_DEVICE_MANAGED amrex::Real cv;
    AMREX_GPU_DEVICE_MANAGED amrex::Real cp;

    AMREX_GPU_DEVICE_MANAGED amrex::Real Pr  = 0.72;     // Prandtl number
    AMREX_GPU_DEVICE_MANAGED amrex::Real C_S = 1.458e-5; // constant in Sutherland's law
    AMREX_GPU_DEVICE_MANAGED amrex::Real T_S = 110.4;    // Sutherland temperature

    AMREX_GPU_DEVICE_MANAGED amrex::Real smallr = 1.e-19;
    AMREX_GPU_DEVICE_MANAGED amrex::Real smallp = 1.e-10;

    void Initialize ()
    {
        constexpr amrex::Real Ru = 8.31451e7;
        cv = Ru / (eos_mu * (eos_gamma-1.0));
        cp = eos_gamma * Ru / (eos_mu * (eos_gamma-1.0));
    }
}
