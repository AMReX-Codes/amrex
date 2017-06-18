
#include <AMReX_Device.H>

bool amrex::Device::in_device_launch_region = false;

#ifdef CUDA
int amrex::Device::cuda_device_id = 0;
#endif

void
amrex::Device::initialize_cuda() {

    initialize_cuda_f();

}
