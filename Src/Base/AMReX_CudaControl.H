#ifndef AMREX_CUDA_CONTROL_H_
#define AMREX_CUDA_CONTROL_H_

namespace amrex {
namespace Cuda {

    extern bool in_launch_region;

    inline bool inLaunchRegion () { return in_launch_region; }

    inline bool setLaunchRegion (bool launch) {
        bool r =  in_launch_region;
        in_launch_region = launch;
        return r;
    }

}
}

#endif
