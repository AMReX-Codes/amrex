#ifndef AMREX_CUDA_DEVICE_H_
#define AMREX_CUDA_DEVICE_H_

#include <cstdlib>
#include <memory>
#include <array>
#include <AMReX.H>
#include <AMReX_Utility.H>
#include <AMReX_GpuError.H>
#include <AMReX_GpuControl.H>

namespace amrex {
namespace Cuda {

class Device
{

public:

    static void Initialize ();
    static void Finalize ();

#if defined(AMREX_USE_CUDA)
    static cudaStream_t cudaStream () noexcept { return cuda_stream; }
#endif
    static void setStreamIndex (const int idx) noexcept;
    static void resetStreamIndex () noexcept { setStreamIndex(-1); }

    static int deviceId () noexcept;

    static void synchronize ();
    static void streamSynchronize ();

    static void htod_memcpy       (void* p_d, const void* p_h, const std::size_t sz);
    static void dtoh_memcpy       (void* p_h, const void* p_d, const std::size_t sz);
    static void htod_memcpy_async (void* p_d, const void* p_h, const std::size_t sz);
    static void dtoh_memcpy_async (void* p_h,const  void* p_d, const std::size_t sz);

    static void mem_advise_set_preferred (void* p, const std::size_t sz, const int device);
    static void mem_advise_set_readonly (void* p, const std::size_t sz);

#if defined(AMREX_USE_CUDA)
    static void setNumThreadsMin (int nx, int ny, int nz) noexcept;
    static void n_threads_and_blocks (const long N, dim3& numBlocks, dim3& numThreads) noexcept;
    static void c_comps_threads_and_blocks (const int* lo, const int* hi, const int comps,
                                            dim3& numBlocks, dim3& numThreads) noexcept;
    static void c_threads_and_blocks (const int* lo, const int* hi, dim3& numBlocks, dim3& numThreads) noexcept;
    static void grid_stride_threads_and_blocks (dim3& numBlocks, dim3& numThreads) noexcept;
    static void box_threads_and_blocks (const Box& bx, dim3& numBlocks, dim3& numThreads) noexcept;

    static std::size_t totalGlobalMem () noexcept { return device_prop.totalGlobalMem; }
    static int numMultiProcessors () noexcept { return device_prop.multiProcessorCount; }
    static int maxThreadsPerMultiProcessor () noexcept { return device_prop.maxThreadsPerMultiProcessor; }
    static int maxThreadsPerBlock () noexcept { return device_prop.maxThreadsPerBlock; }
    static int maxThreadsPerBlock (int dir) noexcept { return device_prop.maxThreadsDim[dir]; }
    static int maxBlocksPerGrid (int dir) noexcept { return device_prop.maxGridSize[dir]; }
#endif

    static std::size_t freeMemAvailable ();

private:

    static void initialize_cuda ();

    static int device_id;
    static int verbose;

#if defined(AMREX_USE_CUDA)
    static constexpr int max_cuda_streams = 16;
    static std::array<cudaStream_t,max_cuda_streams> cuda_streams;
    static cudaStream_t cuda_stream;

    static dim3 numThreadsMin;
    static dim3 numBlocksOverride, numThreadsOverride;

    static cudaDeviceProp device_prop;
#endif
};

}}

#endif
