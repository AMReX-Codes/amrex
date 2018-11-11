
#include <AMReX_GpuControl.H>

namespace amrex {
namespace Cuda {

#if defined(AMREX_USE_GPU)
    bool in_launch_region = true;
#endif

}
}
