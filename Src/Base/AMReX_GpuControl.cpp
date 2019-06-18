
#include <AMReX_GpuControl.H>

namespace amrex {
namespace Gpu {

#if defined(AMREX_USE_GPU)
    bool in_launch_region = true;
    bool in_graph_region = false;
#endif

}
}
