#ifndef AMREX_GPU_ERROR_H_
#define AMREX_GPU_ERROR_H_

#include <AMReX_Utility.H>
#include <string>

#if defined(AMREX_USE_CUDA) && !defined(AMREX_GPU_NO_ERROR_CHECK)

namespace amrex { namespace Cuda {
    inline void ErrorCheck ()
    {
        cudaError_t err = cudaGetLastError();
        if (cudaSuccess != err) {
            std::string errStr(std::string("Cuda last error detected in file ") + __FILE__
                               + " line " + std::to_string(__LINE__)
                               + " " + std::string(cudaGetErrorString(err)));
            amrex::Abort(errStr);
        }
    }
}}

#define AMREX_GPU_SAFE_CALL(call) { \
    cudaError_t amrex_i_err = call; \
    if (cudaSuccess != amrex_i_err) { \
        std::string errStr(std::string("Cuda error in file ") + __FILE__ \
                           + " line " + std::to_string(__LINE__) \
                           + " " + cudaGetErrorString(amrex_i_err)); \
        amrex::Abort(errStr); \
    }}

#define AMREX_GPU_ERROR_CHECK() amrex::Cuda::ErrorCheck()

#else

#define AMREX_GPU_SAFE_CALL(call)  (call)
#define AMREX_GPU_ERROR_CHECK()    ((void)0)

#endif

#endif
