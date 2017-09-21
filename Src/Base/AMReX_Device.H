
#ifndef BL_DEVICE_H
#define BL_DEVICE_H

#include <cstdlib>
#include <memory>

#ifdef NVML
#include <nvml.h>
#endif

extern "C" {
#ifdef AMREX_USE_CUDA
    void initialize_cuda(const int* id, const int* rank);
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
    void gpu_start_profiler();
    void gpu_stop_profiler();
    void set_threads_and_blocks(const int* lo, const int* hi);
    void get_threads_and_blocks(const int* lo, const int* hi, int* bx, int* by, int* bz, int* tx, int* ty, int* tz);
    void gpu_htod_memcpy_async(void* p_d, const void* p_h, const std::size_t* sz, const int* idx);
    void gpu_dtoh_memcpy_async(void* p_h, const void* p_d, const std::size_t* sz, const int* idx);
    void gpu_htod_memprefetch_async(void* p, const std::size_t* sz, const int* idx);
    void gpu_dtoh_memprefetch_async(void* p, const std::size_t* sz, const int* idx);
    void mem_advise_set_preferred(void* p, const std::size_t sz, const int* device);
    void mem_advise_set_readonly(void* p, const std::size_t sz);
    void check_for_gpu_errors();
    void set_is_program_running(int r);
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

    static void stream_synchronize(const int idx = -1);

    static int program_running();

    static void* device_malloc(const std::size_t sz);

    static void device_free(void* ptr);

    static int isProgramRunning();

#ifdef AMREX_USE_CUDA
    template<typename T>
    static std::shared_ptr<T> create_host_pointer(const T* ptr, const int n = 1, const int m = -1);

    template<typename T>
    static std::shared_ptr<T> create_device_pointer(const int n = 1);
#endif

    static void device_htod_memcpy_async(void* p_d, const void* p_h, const std::size_t sz, const int idx);
    static void device_dtoh_memcpy_async(void* p_h, const void* p_d, const std::size_t sz, const int idx);

    static void start_profiler();
    static void stop_profiler();

#if (defined(AMREX_USE_CUDA) && defined(__CUDACC__))
    static void c_threads_and_blocks(const int* lo, const int* hi, dim3& numBlocks, dim3& numThreads);

    static cudaStream_t stream_from_index (int idx);

    static void initialize_cuda_c ();
#endif

private:

    static bool in_device_launch_region;
    static int device_id;
    static int is_device_program_running;

#ifdef AMREX_USE_CUDA
    static const int max_cuda_streams = 100;
    static cudaStream_t cuda_streams[max_cuda_streams];
#endif
};

}

#ifdef AMREX_USE_CUDA
template<typename T>
std::shared_ptr<T>
amrex::Device::create_host_pointer(const T* ptr, const int n, const int m)
{
    std::shared_ptr<T> h_ptr;

    T* temp_ptr;
    const std::size_t sz = n * sizeof(T);
    gpu_hostalloc((void**) &temp_ptr, &sz);

    // Copy only m values if the argument was specified.

    int num;

    if (m < 0)
        num = n;
    else
        num = m;

    for (int i = 0; i < num; ++i)
        *(temp_ptr + i) = *(ptr + i);
    h_ptr.reset(temp_ptr, [](T* ptr) { gpu_freehost(ptr); });

    return h_ptr;
}
#endif



#ifdef AMREX_USE_CUDA
template<typename T>
std::shared_ptr<T>
amrex::Device::create_device_pointer(int n)
{

    std::shared_ptr<T> d_ptr;

    T* temp_ptr = (T*) device_malloc(n * sizeof(T));

    d_ptr.reset(temp_ptr, [](T* ptr) { gpu_free(ptr); });

    return d_ptr;

}
#endif

#endif
