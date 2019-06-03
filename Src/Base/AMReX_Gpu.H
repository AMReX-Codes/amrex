#ifndef AMREX_GPU_H_
#define AMREX_GPU_H_

#ifdef AMREX_USE_HIP
#include <hip/hip_runtime.h>
namespace amrex { namespace Hip {} }
#elif AMREX_USE_CUDA
namespace amrex { namespace Cuda {} }
#endif

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuError.H>

#include <AMReX_GpuUtility.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuAsyncArray.H>
#include <AMReX_GpuAsyncFab.H>

#include <AMReX_CudaMemory.H>
#include <AMReX_CudaRange.H>
#include <AMReX_CudaReduce.H>
#include <AMReX_CudaAllocators.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_CudaElixir.H>
#include <AMReX_CudaGraph.H>

namespace amrex {
    namespace Gpu {
#ifdef AMREX_USE_HIP
        using namespace Hip;
#endif
        using namespace Cuda;
    }
}

#endif
