#ifndef AMREX_GPU_ERROR_H_
#define AMREX_GPU_ERROR_H_

#include <AMReX_Utility.H>
#include <string>

#ifdef AMREX_USE_HIP
#include <hip/hip_runtime.h>
#endif

#if defined(AMREX_USE_GPU) && !defined(AMREX_GPU_NO_ERROR_CHECK)

namespace amrex {

#ifdef AMREX_USE_HIP

using gpuError_t = hipError_t;
constexpr gpuError_t gpuSuccess = hipSuccess;

inline gpuError_t
gpuGetLastError () {
    return hipGetLastError();
}

inline const char*
gpuGetErrorString (gpuError_t error) {
    return hipGetErrorString(error);
}

#else

using gpuError_t = cudaError_t;
constexpr gpuError_t gpuSuccess = cudaSuccess;

inline gpuError_t
gpuGetLastError () {
    return cudaGetLastError();
}

inline const char*
gpuGetErrorString (gpuError_t error) {
    return cudaGetErrorString(error);
}

#endif

namespace Gpu {

    inline void ErrorCheck (const char* file, int line) noexcept
    {
        gpuError_t err = gpuGetLastError();
        if (gpuSuccess != err) {
            std::string errStr(std::string("GPU last error detected in file ") + file
                               + " line " + std::to_string(line)
                               + ": " + std::string(gpuGetErrorString(err)));
            amrex::Abort(errStr);
        }
    }
}
}

#define AMREX_GPU_SAFE_CALL(call) { \
    amrex::gpuError_t amrex_i_err = call;          \
    if (amrex::gpuSuccess != amrex_i_err) { \
        std::string errStr(std::string("GPU error in file ") + __FILE__ \
                           + " line " + std::to_string(__LINE__) \
                           + " " + amrex::gpuGetErrorString(amrex_i_err)); \
        amrex::Abort(errStr); \
    }}

#ifdef AMREX_USE_CUDA
#define AMREX_CUDA_SAFE_CALL(call) { \
    cudaError_t amrex_i_err = call;          \
    if (cudaSuccess != amrex_i_err) { \
        std::string errStr(std::string("CUDA error in file ") + __FILE__ \
                           + " line " + std::to_string(__LINE__) \
                           + " " + cudaGetErrorString(amrex_i_err)); \
        amrex::Abort(errStr); \
    }}
#endif

#ifdef AMREX_USE_HIP
#define AMREX_HIP_SAFE_CALL(call) { \
    hipError_t amrex_i_err = call;          \
    if (hipSuccess != amrex_i_err) { \
        std::string errStr(std::string("HIP error in file ") + __FILE__ \
                           + " line " + std::to_string(__LINE__) \
                           + " " + hipGetErrorString(amrex_i_err)); \
        amrex::Abort(errStr); \
    }}
#endif

#define AMREX_GPU_ERROR_CHECK() amrex::Gpu::ErrorCheck(__FILE__, __LINE__)

#else

#define AMREX_GPU_SAFE_CALL(call)  (call)
#define AMREX_GPU_ERROR_CHECK()    ((void)0)

#endif

#endif
