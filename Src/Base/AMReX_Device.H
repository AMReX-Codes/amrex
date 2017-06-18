
#ifndef BL_DEVICE_H
#define BL_DEVICE_H

#include <cstdlib>

#ifdef CUDA
extern "C" {
    void initialize_cuda();
    void finalize_cuda();
    void gpu_malloc(void** p, const std::size_t* sz);
    void gpu_hostalloc(void** p, const std::size_t* sz);
    void gpu_malloc_managed(void** p, const std::size_t* sz);
    void gpu_free(void* p);
    void gpu_freehost(void* p);
    void gpu_host_device_ptr(void** x, void* y);
    void gpu_synchronize();
    void gpu_htod_memcpy_async(void* p_d, void* p_h, const std::size_t* sz, const int* idx);
    void gpu_dtoh_memcpy_async(void* p_h, void* p_d, const std::size_t* sz, const int* idx);
    void gpu_htod_memprefetch_async(void* p, const std::size_t* sz, const int* idx);
    void gpu_dtoh_memprefetch_async(void* p, const std::size_t* sz, const int* idx);
    void mem_advise_set_preferred(void* p, const std::size_t sz, const int* device);
    void mem_advise_set_readonly(void* p, const std::size_t sz);
}
#endif

namespace amrex {

class Device
{

public:

    static void beginDeviceLaunchRegion() { in_device_launch_region = true; }

    static void endDeviceLaunchRegion() { in_device_launch_region = false; }

    static void setDeviceLaunchRegion(bool r) { in_device_launch_region = r; }

    static bool inDeviceLaunchRegion() {
#ifdef CUDA
      gpu_synchronize();
#endif
      return in_device_launch_region; 
    }

#ifdef CUDA
    static int cudaDeviceId() { return cuda_device_id; }
#endif

private:

    static bool in_device_launch_region;
    static int cuda_device_id;

};

}

#endif
