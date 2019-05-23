#ifndef AMREX_GPU_H_
#define AMREX_GPU_H_

#ifdef AMREX_USE_HIP
#include <hip/hip_runtime.h>
#endif

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuError.H>

#include <AMReX_CudaMemory.H>
#include <AMReX_CudaRange.H>
#include <AMReX_CudaUtility.H>
#include <AMReX_CudaReduce.H>
#include <AMReX_CudaDevice.H>
#include <AMReX_CudaAllocators.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_CudaAsyncFab.H>
#include <AMReX_CudaAsyncArray.H>
#include <AMReX_CudaElixir.H>

namespace amrex {
    namespace Gpu {
        using namespace Cuda;
    }
}

#endif
