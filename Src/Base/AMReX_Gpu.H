#ifndef AMREX_GPU_H_
#define AMREX_GPU_H_
#include <AMReX_Config.H>

#ifdef AMREX_USE_HIP
#include <hip/hip_runtime.h>
namespace amrex::Hip {}
#elif defined(AMREX_USE_CUDA)
namespace amrex::Cuda {}
#endif

#include <AMReX_GpuQualifiers.H> // IWYU pragma: export
#include <AMReX_GpuKernelInfo.H> // IWYU pragma: export
#include <AMReX_GpuPrint.H> // IWYU pragma: export
#include <AMReX_GpuAssert.H> // IWYU pragma: export
#include <AMReX_GpuTypes.H> // IWYU pragma: export
#include <AMReX_GpuControl.H> // IWYU pragma: export
#include <AMReX_GpuLaunch.H> // IWYU pragma: export
#include <AMReX_GpuError.H> // IWYU pragma: export

#include <AMReX_GpuAtomic.H> // IWYU pragma: export
#include <AMReX_GpuUtility.H> // IWYU pragma: export
#include <AMReX_GpuDevice.H> // IWYU pragma: export
#include <AMReX_GpuBuffer.H> // IWYU pragma: export
#include <AMReX_GpuAsyncArray.H> // IWYU pragma: export
#include <AMReX_GpuElixir.H> // IWYU pragma: export
#include <AMReX_GpuMemory.H> // IWYU pragma: export
#include <AMReX_GpuRange.H> // IWYU pragma: export
#include <AMReX_GpuReduce.H> // IWYU pragma: export

#include <AMReX_GpuAllocators.H> // IWYU pragma: export
#include <AMReX_GpuContainers.H> // IWYU pragma: export
#include <AMReX_CudaGraph.H> // IWYU pragma: export

namespace amrex::Gpu {
#ifdef AMREX_USE_HIP
        using namespace Hip;
#elif defined(AMREX_USE_CUDA)
        using namespace Cuda;
#endif
}

#endif
