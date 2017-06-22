
#include <AMReX_Device.H>

bool amrex::Device::in_device_launch_region = false;

#ifdef CUDA
int amrex::Device::cuda_device_id = 0;
#endif

void
amrex::Device::initialize_cuda() {

#ifdef CUDA
    initialize_cuda_f();

    get_cuda_device_id(&cuda_device_id);
#endif

}

void
amrex::Device::set_stream_index(const int idx) {

#ifdef CUDA
    set_stream_idx(idx);
#endif

}

int
amrex::Device::get_stream_index() {

    int index = -1;

#ifdef CUDA
    get_stream_idx(&index);
#endif

    return index;

}

void*
amrex::Device::get_host_pointer(const void* ptr) {

    void* r;
#ifdef CUDA
    gpu_host_device_ptr(&r, ptr);
#endif
    return r;

}

void
amrex::Device::check_for_errors() {

#ifdef CUDA
    check_for_gpu_errors();
#endif

}
