
#include <AMReX_GpuControl.H>

namespace amrex::Gpu {

#if defined(AMREX_USE_GPU)
bool in_launch_region = true;
bool in_graph_region = false;
bool in_single_stream_region = false;
bool in_nosync_region = false;
#endif

}
