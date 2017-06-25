
#ifndef BL_DEVICE_H
#define BL_DEVICE_H

#include <cstdlib>

extern "C" {
#ifdef CUDA
    void initialize_cuda();
    void finalize_cuda();
    void set_stream_idx(const int idx);
    void get_stream_idx(int* idx);
    void get_cuda_device_id(int* device);
    void gpu_malloc(void** p, const std::size_t* sz);
    void gpu_hostalloc(void** p, const std::size_t* sz);
    void gpu_malloc_managed(void** p, const std::size_t* sz);
    void gpu_free(void* p);
    void gpu_freehost(void* p);
    void gpu_host_device_ptr(void** x, const void* y);
    void gpu_stream_synchronize(const int idx);
    void gpu_synchronize();
    void set_threads_and_blocks(const int* lo, const int* hi);
    void gpu_htod_memcpy_async(void* p_d, void* p_h, const std::size_t* sz, const int* idx);
    void gpu_dtoh_memcpy_async(void* p_h, void* p_d, const std::size_t* sz, const int* idx);
    void gpu_htod_memprefetch_async(void* p, const std::size_t* sz, const int* idx);
    void gpu_dtoh_memprefetch_async(void* p, const std::size_t* sz, const int* idx);
    void mem_advise_set_preferred(void* p, const std::size_t sz, const int* device);
    void mem_advise_set_readonly(void* p, const std::size_t sz);
    void check_for_gpu_errors();
#endif
}

namespace amrex {

class Device
{

public:

    static void initialize_device();

    static void finalize_device();

    static void set_stream_index(const int idx);

    static int  get_stream_index();

    static void prepare_for_launch(const int* lo, const int* hi);

    static void beginDeviceLaunchRegion() { in_device_launch_region = true; }

    static void endDeviceLaunchRegion() { in_device_launch_region = false; }

    static void setDeviceLaunchRegion(bool r) { in_device_launch_region = r; }

    static bool inDeviceLaunchRegion() { return in_device_launch_region; }

    static void* get_host_pointer(const void* ptr);

    static int deviceId();

    static void check_for_errors();

    static void synchronize();

    static void* device_malloc(const std::size_t sz);

    static void device_free(void* ptr);

    static void device_htod_memcpy_async(void* p_d, void* p_h, const std::size_t sz, const int idx);
    static void device_dtoh_memcpy_async(void* p_h, void* p_d, const std::size_t sz, const int idx);

private:

    static bool in_device_launch_region;
    static int device_id;

};

}

#endif
