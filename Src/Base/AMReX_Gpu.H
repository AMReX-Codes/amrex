#ifndef AMREX_GPU_H_
#define AMREX_GPU_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>

#include <AMReX_CudaMemory.H>
#include <AMReX_CudaUtility.H>
#include <AMReX_CudaDevice.H>
#include <AMReX_CudaAllocators.H>
#include <AMReX_CudaContainers.H>

namespace amrex {
    namespace Gpu = Cuda;
}

#endif
