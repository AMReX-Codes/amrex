
#include <CNS_parm.H>

namespace Parm
{
    AMREX_GPU_DEVICE_MANAGED amrex::Real eos_gamma = 1.4;

    AMREX_GPU_DEVICE_MANAGED amrex::Real smallr = 1.e-19;
    AMREX_GPU_DEVICE_MANAGED amrex::Real smallp = 1.e-10;
}
