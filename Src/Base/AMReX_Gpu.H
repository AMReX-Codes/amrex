#ifndef AMREX_GPU_H_
#define AMREX_GPU_H_

#ifdef AMREX_USE_HIP
#include <hip/hip_runtime.h>
namespace amrex { namespace Hip {} }
#elif defined(AMREX_USE_CUDA)
namespace amrex { namespace Cuda {} }
#endif

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuPrint.H>
#include <AMReX_GpuAssert.H>
#include <AMReX_GpuTypes.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuFuse.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuError.H>

#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuUtility.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuAsyncArray.H>
#include <AMReX_GpuElixir.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_GpuRange.H>
#include <AMReX_GpuReduce.H>

#include <AMReX_GpuAllocators.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_CudaGraph.H>

namespace amrex {
    namespace Gpu {
#ifdef AMREX_USE_HIP
        using namespace Hip;
#elif defined(AMREX_USE_CUDA)
        using namespace Cuda;
#endif
    }
}

#endif
