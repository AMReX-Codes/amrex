
#include <AMReX_GpuControl.H>

namespace amrex {
namespace Gpu {

#if defined(AMREX_USE_GPU)
    bool in_launch_region = true;
#else
    bool in_launch_region = false;
#endif

}
}
