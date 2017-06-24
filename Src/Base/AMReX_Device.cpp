
#include <AMReX_Device.H>

bool amrex::Device::in_device_launch_region = false;
int amrex::Device::device_id = 0;

void
amrex::Device::initialize_device() {

#ifdef CUDA
    initialize_cuda();

    get_cuda_device_id(&device_id);
#endif

}

void
amrex::Device::finalize_device() {

#ifdef CUDA
    finalize_cuda();
#endif

}

int
amrex::Device::deviceId() {

    return device_id;

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

void
amrex::Device::prepare_for_launch(const int* lo, const int* hi) {

    // Sets the number of threads and blocks in Fortran.

#ifdef CUDA
    set_threads_and_blocks(lo, hi);
#endif

}

void*
amrex::Device::get_host_pointer(const void* ptr) {

    void* r = NULL;
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

void
amrex::Device::synchronize() {

#ifdef CUDA
    gpu_synchronize();
#endif

}
