
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

    void* r = const_cast<void*>(ptr);
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

void
amrex::Device::stream_synchronize(const int idx) {

#ifdef CUDA
    gpu_stream_synchronize(idx);
#endif

}

void*
amrex::Device::device_malloc(const std::size_t sz) {

    void* ptr = NULL;

#ifdef CUDA
    gpu_malloc(&ptr, &sz);
#endif

    return ptr;

}

void
amrex::Device::device_free(void* ptr) {

#ifdef CUDA
    gpu_free(ptr);
#endif

}

void
amrex::Device::device_htod_memcpy_async(void* p_d, const void* p_h, const std::size_t sz, const int idx) {

#ifdef CUDA
    gpu_htod_memcpy_async(p_d, p_h, &sz, &idx);
#endif

}

void
amrex::Device::device_dtoh_memcpy_async(void* p_h, const void* p_d, const std::size_t sz, const int idx) {

#ifdef CUDA
    gpu_dtoh_memcpy_async(p_h, p_d, &sz, &idx);
#endif

}

void
amrex::Device::start_profiler() {

#ifdef CUDA
    gpu_start_profiler();
#endif

}

void
amrex::Device::stop_profiler() {

#ifdef CUDA
    gpu_stop_profiler();
#endif

}
