#ifndef AMREX_GPU_KERNEL_INFO_H_
#define AMREX_GPU_KERNEL_INFO_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_GpuKernelInfo.H
 * \brief Metadata tag passed to GPU kernel launchers to describe kernel properties.
 */

namespace amrex::Gpu {

/**
 * \class amrex::Gpu::KernelInfo
 * \ingroup amrex_execution
 * \brief Metadata tag describing optional properties of a GPU kernel launch.
 *
 * Pass a \c KernelInfo object to kernel-launch helpers (e.g., \c ParallelFor)
 * to advertise kernel properties such as whether the kernel performs a
 * device-side reduction.
 */
class KernelInfo
{
public:
    //! Set whether the kernel performs a device-side reduction; returns \c *this for chaining.
    KernelInfo& setReduction (bool flag) { has_reduction = flag; return *this; }
    //! Return \c true if a device-side reduction was requested.
    [[nodiscard]] bool hasReduction () const { return has_reduction; }
private:
    bool has_reduction = false;
};

}

#endif
