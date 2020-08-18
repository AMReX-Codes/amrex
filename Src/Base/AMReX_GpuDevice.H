#ifndef AMREX_GPU_DEVICE_H_
#define AMREX_GPU_DEVICE_H_

#include <cstdlib>
#include <memory>
#include <array>
#include <AMReX.H>
#include <AMReX_Utility.H>
#include <AMReX_GpuTypes.H>
#include <AMReX_GpuError.H>
#include <AMReX_GpuControl.H>
#include <AMReX_OpenMP.H>
#include <AMReX_Vector.H>

namespace amrex {

#ifdef AMREX_USE_HIP
using gpuDeviceProp_t = hipDeviceProp_t;
#elif defined(AMREX_USE_CUDA)
using gpuDeviceProp_t = cudaDeviceProp;
#elif defined(AMREX_USE_DPCPP)
    struct gpuDeviceProp_t {
        std::string name;
        std::size_t totalGlobalMem;
        std::size_t sharedMemPerBlock;
        int multiProcessorCount;
        int maxThreadsPerMultiProcessor;
        int maxThreadsPerBlock;
        int maxThreadsDim[3];
        int maxGridSize[3];
        int warpSize;
        Long maxMemAllocSize; // oneAPI only
        int managedMemory;
        int concurrentManagedAccess;
        int maxParameterSize;
    };
#endif

namespace Gpu {

class Device
{

public:

    static void Initialize ();
    static void Finalize ();

#if defined(AMREX_USE_GPU)
    static gpuStream_t gpuStream () noexcept { return gpu_stream[OpenMP::get_thread_num()]; }
    static gpuStream_t nullStream () noexcept { return gpu_default_stream; }
#ifdef AMREX_USE_CUDA
    // for backward compatibility
    static cudaStream_t cudaStream () noexcept { return gpu_stream[OpenMP::get_thread_num()]; }
#endif
#ifdef AMREX_USE_DPCPP
    static sycl::queue& nullQueue () noexcept { return *(gpu_default_stream.queue); }
    static sycl::queue& streamQueue () noexcept { return *(gpu_stream[OpenMP::get_thread_num()].queue); }
    static sycl::queue& streamQueue (int i) noexcept { return *(gpu_streams[i].queue); }
    static bool onNullStream () noexcept { return gpu_stream[OpenMP::get_thread_num()] == gpu_default_stream; }
    static bool onNullStream (gpuStream_t stream) noexcept { return stream == gpu_default_stream; }
#endif
#endif

    static int numGpuStreams () noexcept { return max_gpu_streams; }
    static void setStreamIndex (const int idx) noexcept;
    static void resetStreamIndex () noexcept { setStreamIndex(-1); }

#ifdef AMREX_USE_GPU
    static gpuStream_t setStream (gpuStream_t s) noexcept;
    static gpuStream_t resetStream () noexcept;
#endif

    static int deviceId () noexcept;
    static int numDevicesUsed () noexcept;

    static void synchronize () noexcept;
    static void streamSynchronize () noexcept;
#ifdef AMREX_USE_DPCPP
    static void nonNullStreamSynchronize () noexcept;
#endif

#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10) )
    // Generic graph selection. These should be called by users.
    static void startGraphRecording(bool first_iter, void* h_ptr, void* d_ptr, size_t sz);
    static cudaGraphExec_t stopGraphRecording(bool last_iter);

    // Instantiate a created cudaGtaph
    static cudaGraphExec_t instantiateGraph(cudaGraph_t graph);

    // Execute an instantiated cudaGraphExec
    static void executeGraph(const cudaGraphExec_t &graphExec, bool synch = true);

#endif

    static void mem_advise_set_preferred (void* p, const std::size_t sz, const int device);
    static void mem_advise_set_readonly (void* p, const std::size_t sz);

#ifdef AMREX_USE_GPU
    static void setNumThreadsMin (int nx, int ny, int nz) noexcept;
    static void n_threads_and_blocks (const Long N, dim3& numBlocks, dim3& numThreads) noexcept;
    static void c_comps_threads_and_blocks (const int* lo, const int* hi, const int comps,
                                            dim3& numBlocks, dim3& numThreads) noexcept;
    static void c_threads_and_blocks (const int* lo, const int* hi, dim3& numBlocks, dim3& numThreads) noexcept;
    static void grid_stride_threads_and_blocks (dim3& numBlocks, dim3& numThreads) noexcept;
    static void box_threads_and_blocks (const Box& bx, dim3& numBlocks, dim3& numThreads) noexcept;

    static std::size_t totalGlobalMem () noexcept { return device_prop.totalGlobalMem; }
    static std::size_t sharedMemPerBlock () noexcept { return device_prop.sharedMemPerBlock; }
    static int numMultiProcessors () noexcept { return device_prop.multiProcessorCount; }
    static int maxThreadsPerMultiProcessor () noexcept { return device_prop.maxThreadsPerMultiProcessor; }
    static int maxThreadsPerBlock () noexcept { return device_prop.maxThreadsPerBlock; }
    static int maxThreadsPerBlock (int dir) noexcept { return device_prop.maxThreadsDim[dir]; }
    static int maxBlocksPerGrid (int dir) noexcept { return device_prop.maxGridSize[dir]; }
    static std::string deviceName () noexcept { return std::string(device_prop.name); }
#endif

    static std::size_t freeMemAvailable ();

#ifdef AMREX_USE_GPU
    static constexpr int warp_size = AMREX_HIP_OR_CUDA_OR_DPCPP(64,32,16);

    static int maxBlocksPerLaunch () noexcept { return max_blocks_per_launch; }

#ifdef AMREX_USE_DPCPP
    static Long maxMemAllocSize () noexcept { return device_prop.maxMemAllocSize; }
    static sycl::context& syclContext () { return *sycl_context; }
    static sycl::device& syclDevice () { return *sycl_device; }
#endif
#endif

private:

    static void initialize_gpu ();

    static int device_id;
    static int num_devices_used;
    static int verbose;

#ifdef AMREX_USE_GPU
    static constexpr int max_gpu_streams = 4;
#else
    // Equivalent to "single dependent stream". Fits best
    //  with math this is used in ("x/max_streams").
    static constexpr int max_gpu_streams = 1;
#endif

#ifdef AMREX_USE_GPU
    static dim3 numThreadsMin;
    static dim3 numBlocksOverride, numThreadsOverride;

    static std::array<gpuStream_t,max_gpu_streams> gpu_streams;
    static gpuStream_t gpu_default_stream;
    static Vector<gpuStream_t> gpu_stream;
    static gpuDeviceProp_t device_prop;
    static int max_blocks_per_launch;

#ifdef AMREX_USE_DPCPP
    static std::unique_ptr<sycl::context> sycl_context;
    static std::unique_ptr<sycl::device>  sycl_device;
#endif
#endif
};

// Put these in amrex::Gpu

#if defined(AMREX_USE_GPU)
inline gpuStream_t
gpuStream () noexcept
{
    return Device::gpuStream();
}

inline gpuStream_t
nullStream () noexcept
{
    return Device::nullStream();
}
#endif

inline int
numGpuStreams () noexcept
{
    return Device::numGpuStreams();
}

inline void
synchronize () noexcept
{
    Device::synchronize();
}

inline void
streamSynchronize () noexcept
{
    Device::streamSynchronize();
}

#ifdef AMREX_USE_DPCPP
inline void
nonNullStreamSynchronize () noexcept
{
    Device::nonNullStreamSynchronize();
}
#endif

#ifdef AMREX_USE_GPU

inline void
htod_memcpy (void* p_d, const void* p_h, const std::size_t sz) noexcept
{
    if (sz == 0) return;
#ifdef AMREX_USE_DPCPP
    Device::nonNullStreamSynchronize();
    auto& q = Device::nullQueue();
    q.submit([&] (sycl::handler& h) { h.memcpy(p_d, p_h, sz); });
    try {
        q.wait_and_throw();
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("htod_memcpy: ")+ex.what()+"!!!!!");
    }
#else
    AMREX_HIP_OR_CUDA(
        AMREX_HIP_SAFE_CALL(hipMemcpy(p_d, p_h, sz, hipMemcpyHostToDevice));,
        AMREX_CUDA_SAFE_CALL(cudaMemcpy(p_d, p_h, sz, cudaMemcpyHostToDevice)); )
#endif
}

inline void
dtoh_memcpy (void* p_h, const void* p_d, const std::size_t sz) noexcept
{
    if (sz == 0) return;
#ifdef AMREX_USE_DPCPP
    Device::nonNullStreamSynchronize();
    auto& q = Device::nullQueue();
    q.submit([&] (sycl::handler& h) { h.memcpy(p_h, p_d, sz); });
    try {
        q.wait_and_throw();
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("dtoh_memcpy: ")+ex.what()+"!!!!!");
    }
    Gpu::synchronize(); // To mimic cuda behavior
#else
    AMREX_HIP_OR_CUDA(
        AMREX_HIP_SAFE_CALL(hipMemcpy(p_h, p_d, sz, hipMemcpyDeviceToHost));,
        AMREX_CUDA_SAFE_CALL(cudaMemcpy(p_h, p_d, sz, cudaMemcpyDeviceToHost)); )
#endif
}

inline void
dtod_memcpy (void* p_d_dst, const void* p_d_src, const std::size_t sz) noexcept
{
    if (sz == 0) return;
#ifdef AMREX_USE_DPCPP
    Device::nonNullStreamSynchronize();
    auto& q = Device::nullQueue();
    q.submit([&] (sycl::handler& h) { h.memcpy(p_d_dst, p_d_src, sz); });
    try {
        q.wait_and_throw();
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("dtod_memcpy: ")+ex.what()+"!!!!!");
    }
#else
    AMREX_HIP_OR_CUDA(
        AMREX_HIP_SAFE_CALL(hipMemcpy(p_d_dst, p_d_src, sz, hipMemcpyDeviceToDevice));,
        AMREX_CUDA_SAFE_CALL(cudaMemcpy(p_d_dst, p_d_src, sz, cudaMemcpyDeviceToDevice)); )
#endif
}

inline void
htod_memcpy_async (void* p_d, const void* p_h, const std::size_t sz) noexcept
{
    if (sz == 0) return;
#ifdef AMREX_USE_DPCPP
    auto& q = Device::streamQueue();
    q.submit([&] (sycl::handler& h) { h.memcpy(p_d, p_h, sz); });
#else
    AMREX_HIP_OR_CUDA(
        AMREX_HIP_SAFE_CALL(hipMemcpyAsync(p_d, p_h, sz, hipMemcpyHostToDevice, gpuStream()));,
        AMREX_CUDA_SAFE_CALL(cudaMemcpyAsync(p_d, p_h, sz, cudaMemcpyHostToDevice, gpuStream())); )
#endif
}

inline void
dtoh_memcpy_async (void* p_h, const void* p_d, const std::size_t sz) noexcept
{
    if (sz == 0) return;
#ifdef AMREX_USE_DPCPP
    auto& q = Device::streamQueue();
    q.submit([&] (sycl::handler& h) { h.memcpy(p_h, p_d, sz); });
#else
    AMREX_HIP_OR_CUDA(
        AMREX_HIP_SAFE_CALL(hipMemcpyAsync(p_h, p_d, sz, hipMemcpyDeviceToHost, gpuStream()));,
        AMREX_CUDA_SAFE_CALL(cudaMemcpyAsync(p_h, p_d, sz, cudaMemcpyDeviceToHost, gpuStream())); )
#endif
}

inline void
dtod_memcpy_async (void* p_d_dst, const void* p_d_src, const std::size_t sz) noexcept
{
    if (sz == 0) return;
#ifdef AMREX_USE_DPCPP
    auto& q = Device::streamQueue();
    q.submit([&] (sycl::handler& h) { h.memcpy(p_d_dst, p_d_src, sz); });
#else
    AMREX_HIP_OR_CUDA(
        AMREX_HIP_SAFE_CALL(hipMemcpyAsync(p_d_dst, p_d_src, sz, hipMemcpyDeviceToDevice, gpuStream()));,
        AMREX_CUDA_SAFE_CALL(cudaMemcpyAsync(p_d_dst, p_d_src, sz, cudaMemcpyDeviceToDevice, gpuStream())); )
#endif
}

#endif

#ifdef AMREX_USE_DPCPP
inline bool
onNullStream ()
{
    return Device::onNullStream();
}

inline bool
onNullStream (gpuStream_t stream)
{
    return Device::onNullStream(stream);
}
#endif

#ifdef AMREX_USE_GPU
void callbackAdded ();
void resetNumCallbacks ();
int getNumCallbacks ();
#endif

}}

#endif
