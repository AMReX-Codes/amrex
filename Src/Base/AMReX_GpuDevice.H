#ifndef AMREX_GPU_DEVICE_H_
#define AMREX_GPU_DEVICE_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_Extension.H>
#include <AMReX_Utility.H>
#include <AMReX_GpuTypes.H>
#include <AMReX_GpuError.H>
#include <AMReX_GpuControl.H>
#include <AMReX_OpenMP.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstring>
#include <memory>

#define AMREX_GPU_MAX_STREAMS 8

#ifdef AMREX_USE_GPU
namespace amrex {
#ifdef AMREX_USE_HIP
using gpuDeviceProp_t = hipDeviceProp_t;
#elif defined(AMREX_USE_CUDA)
using gpuDeviceProp_t = cudaDeviceProp;
#elif defined(AMREX_USE_SYCL)
    struct gpuDeviceProp_t {
        std::string name;
        std::string vendor; // SYCL only (inferred for CUDA and HIP)
        std::size_t totalGlobalMem;
        std::size_t sharedMemPerBlock;
        int multiProcessorCount;
        int maxThreadsPerMultiProcessor;
        int maxThreadsPerBlock;
        int maxThreadsDim[3];
        int maxGridSize[3];
        int warpSize;
        Long maxMemAllocSize; // SYCL only
        int managedMemory;
        int concurrentManagedAccess;
        int maxParameterSize;
    };
#endif
}
#endif

namespace amrex::Gpu {

class Device
{

public:

    static void Initialize ();
    static void Finalize ();

#if defined(AMREX_USE_GPU)
    static gpuStream_t gpuStream () noexcept { return gpu_stream[OpenMP::get_thread_num()]; }
#ifdef AMREX_USE_CUDA
    /** for backward compatibility */
    static cudaStream_t cudaStream () noexcept { return gpu_stream[OpenMP::get_thread_num()]; }
#endif
#ifdef AMREX_USE_SYCL
    static sycl::queue& streamQueue () noexcept { return *(gpu_stream[OpenMP::get_thread_num()].queue); }
    static sycl::queue& streamQueue (int i) noexcept { return *(gpu_stream_pool[i].queue); }
#endif
#endif

    static int numGpuStreams () noexcept {
        return inSingleStreamRegion() ? 1 : max_gpu_streams;
    }

    static void setStreamIndex (int idx) noexcept;
    static void resetStreamIndex () noexcept { setStreamIndex(0); }

#ifdef AMREX_USE_GPU
    static int streamIndex (gpuStream_t s = gpuStream()) noexcept;

    static gpuStream_t setStream (gpuStream_t s) noexcept;
    static gpuStream_t resetStream () noexcept;
#endif

    static int deviceId () noexcept;
    static int numDevicesUsed () noexcept; // Total number of device used
    static int numDevicePartners () noexcept; // Number of partners sharing my device

    /**
     * Halt execution of code until GPU has finished processing all previously requested
     * tasks.
     */
    static void synchronize () noexcept;

    /**
     * Halt execution of code until the current AMReX GPU stream has finished processing all
     * previously requested tasks.
     */
    static void streamSynchronize () noexcept;

    /**
     * Halt execution of code until all AMReX GPU streams have finished processing all
     * previously requested tasks.
     */
    static void streamSynchronizeAll () noexcept;

#if defined(__CUDACC__)
    /**  Generic graph selection. These should be called by users.  */
    static void startGraphRecording(bool first_iter, void* h_ptr, void* d_ptr, size_t sz);
    static cudaGraphExec_t stopGraphRecording(bool last_iter);

    /** Instantiate a created cudaGtaph */
    static cudaGraphExec_t instantiateGraph(cudaGraph_t graph);

    /** Execute an instantiated cudaGraphExec */
    static void executeGraph(const cudaGraphExec_t &graphExec, bool synch = true);

#endif

    static void mem_advise_set_preferred (void* p, std::size_t sz, int device);
    static void mem_advise_set_readonly (void* p, std::size_t sz);

#ifdef AMREX_USE_GPU
    static void setNumThreadsMin (int nx, int ny, int nz) noexcept;
    static void n_threads_and_blocks (const Long N, dim3& numBlocks, dim3& numThreads) noexcept;
    static void c_comps_threads_and_blocks (const int* lo, const int* hi, const int comps,
                                            dim3& numBlocks, dim3& numThreads) noexcept;
    static void c_threads_and_blocks (const int* lo, const int* hi, dim3& numBlocks, dim3& numThreads) noexcept;
    static void grid_stride_threads_and_blocks (dim3& numBlocks, dim3& numThreads) noexcept;

    static std::size_t totalGlobalMem () noexcept { return device_prop.totalGlobalMem; }
    static std::size_t sharedMemPerBlock () noexcept { return device_prop.sharedMemPerBlock; }
    static int numMultiProcessors () noexcept { return device_prop.multiProcessorCount; }
    static int maxThreadsPerMultiProcessor () noexcept { return device_prop.maxThreadsPerMultiProcessor; }
    static int maxThreadsPerBlock () noexcept { return device_prop.maxThreadsPerBlock; }
    static int maxThreadsPerBlock (int dir) noexcept { return device_prop.maxThreadsDim[dir]; }
    static int maxBlocksPerGrid (int dir) noexcept { return device_prop.maxGridSize[dir]; }
    static std::string deviceName () noexcept { return std::string(device_prop.name); }
#endif

#ifdef AMREX_USE_CUDA
    static int devicePropMajor () noexcept { return device_prop.major; }
    static int devicePropMinor () noexcept { return device_prop.minor; }
#endif

    static std::string deviceVendor() noexcept
    {
#if defined(AMREX_USE_HIP) && defined(__HIP_PLATFORM_AMD__)
        return std::string("AMD");
#elif defined(AMREX_USE_CUDA) || (defined(AMREX_USE_HIP) && defined(__HIP_PLATFORM_NVIDIA__))
        // Using HIP on NVIDIA GPUs isn't currently supported by AMReX
        return std::string("NVIDIA");
#elif defined(AMREX_USE_SYCL)
        return device_prop.vendor;
#else
        return std::string("Unknown");
#endif
    }

    static std::size_t freeMemAvailable ();
    static void profilerStart ();
    static void profilerStop ();

#ifdef AMREX_USE_GPU

    static int memoryPoolsSupported () noexcept { return memory_pools_supported; }

#if defined(AMREX_USE_HIP)
    static AMREX_EXPORT constexpr int warp_size = AMREX_AMDGCN_WAVEFRONT_SIZE;
#elif defined(AMREX_USE_SYCL)
    static AMREX_EXPORT constexpr int warp_size = AMREX_SYCL_SUB_GROUP_SIZE;
#else
    static AMREX_EXPORT constexpr int warp_size = AMREX_HIP_OR_CUDA(64,32);
#endif

    static unsigned int maxBlocksPerLaunch () noexcept { return max_blocks_per_launch; }

#ifdef AMREX_USE_SYCL
    static Long maxMemAllocSize () noexcept { return device_prop.maxMemAllocSize; }
    static sycl::context& syclContext () { return *sycl_context; }
    static sycl::device& syclDevice () { return *sycl_device; }
#endif
#endif

private:

    static void initialize_gpu ();

    static AMREX_EXPORT int device_id;
    static AMREX_EXPORT int num_devices_used;
    static AMREX_EXPORT int num_device_partners;
    static AMREX_EXPORT int verbose;
    static AMREX_EXPORT int max_gpu_streams;

#ifdef AMREX_USE_GPU
    static AMREX_EXPORT dim3 numThreadsMin;
    static AMREX_EXPORT dim3 numBlocksOverride, numThreadsOverride;

    static AMREX_EXPORT Vector<gpuStream_t> gpu_stream_pool; // The size of this is max_gpu_stream
    // The non-owning gpu_stream is used to store the current stream that will be used.
    // gpu_stream is a vector so that it's thread safe to write to it.
    static AMREX_EXPORT Vector<gpuStream_t> gpu_stream; // The size of this is omp_max_threads
    static AMREX_EXPORT gpuDeviceProp_t device_prop;
    static AMREX_EXPORT int memory_pools_supported;
    static AMREX_EXPORT unsigned int max_blocks_per_launch;

#ifdef AMREX_USE_SYCL
    static AMREX_EXPORT std::unique_ptr<sycl::context> sycl_context;
    static AMREX_EXPORT std::unique_ptr<sycl::device>  sycl_device;
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

inline void
streamSynchronizeAll () noexcept
{
    Device::streamSynchronizeAll();
}

#ifdef AMREX_USE_GPU

inline void
htod_memcpy_async (void* p_d, const void* p_h, const std::size_t sz) noexcept
{
    if (sz == 0) { return; }
#ifdef AMREX_USE_SYCL
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
    if (sz == 0) { return; }
#ifdef AMREX_USE_SYCL
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
    if (sz == 0) { return; }
#ifdef AMREX_USE_SYCL
    auto& q = Device::streamQueue();
    q.submit([&] (sycl::handler& h) { h.memcpy(p_d_dst, p_d_src, sz); });
#else
    AMREX_HIP_OR_CUDA(
        AMREX_HIP_SAFE_CALL(hipMemcpyAsync(p_d_dst, p_d_src, sz, hipMemcpyDeviceToDevice, gpuStream()));,
        AMREX_CUDA_SAFE_CALL(cudaMemcpyAsync(p_d_dst, p_d_src, sz, cudaMemcpyDeviceToDevice, gpuStream())); )
#endif
}

inline void
htod_memcpy (void* p_d, const void* p_h, const std::size_t sz) noexcept
{
    if (sz == 0) { return; }
    htod_memcpy_async(p_d, p_h, sz);
    Gpu::streamSynchronize();
}

inline void
dtoh_memcpy (void* p_h, const void* p_d, const std::size_t sz) noexcept
{
    if (sz == 0) { return; }
    dtoh_memcpy_async(p_h, p_d, sz);
    Gpu::streamSynchronize();
}

inline void
dtod_memcpy (void* p_d_dst, const void* p_d_src, const std::size_t sz) noexcept
{
    if (sz == 0) { return; }
    dtod_memcpy_async(p_d_dst, p_d_src, sz);
    Gpu::streamSynchronize();
}

#endif

#ifdef AMREX_USE_HYPRE
void hypreSynchronize ();
#endif

//! Copy `nbytes` bytes from host to device global variable. `offset` is the
//! offset in bytes from the start of the device global variable.
template <typename T>
void memcpy_from_host_to_device_global_async (T& dg, const void* src,
                                              std::size_t nbytes,
                                              std::size_t offset = 0)
{
#if defined(AMREX_USE_CUDA)
    AMREX_CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(dg, src, nbytes, offset,
                                                 cudaMemcpyHostToDevice,
                                                 Device::gpuStream()));
#elif defined(AMREX_USE_HIP)
    AMREX_HIP_SAFE_CALL(hipMemcpyToSymbolAsync(dg, src, nbytes, offset,
                                               hipMemcpyHostToDevice,
                                               Device::gpuStream()));
#elif defined(AMREX_USE_SYCL)
    Device::streamQueue().memcpy(dg, src, nbytes, offset);
#else
    auto* p = (char*)(&dg);
    std::memcpy(p+offset, src, nbytes);
#endif
}

//! Copy `nbytes` bytes from device global variable to host. `offset` is the
//! offset in bytes from the start of the device global variable.
template <typename T>
void memcpy_from_device_global_to_host_async (void* dst, T const& dg,
                                              std::size_t nbytes,
                                              std::size_t offset = 0)
{
#if defined(AMREX_USE_CUDA)
    AMREX_CUDA_SAFE_CALL(cudaMemcpyFromSymbolAsync(dst, dg, nbytes, offset,
                                                   cudaMemcpyDeviceToHost,
                                                   Device::gpuStream()));
#elif defined(AMREX_USE_HIP)
    AMREX_HIP_SAFE_CALL(hipMemcpyFromSymbolAsync(dst, dg, nbytes, offset,
                                                  hipMemcpyDeviceToHost,
                                                  Device::gpuStream()));
#elif defined(AMREX_USE_SYCL)
    Device::streamQueue().memcpy(dst, dg, nbytes, offset);
#else
    auto const* p = (char const*)(&dg);
    std::memcpy(dst, p+offset, nbytes);
#endif
}

}

#endif
