
#include <AMReX_CudaControl.H>

namespace amrex {
namespace Cuda {

#if defined(AMREX_USE_CUDA)
    bool in_launch_region = true;
#else
    bool in_launch_region = false;
#endif

}
}
