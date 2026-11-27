#ifndef AMREX_GPU_SPARSE_H_
#define AMREX_GPU_SPARSE_H_
#include <AMReX_Config.H>

#if defined(AMREX_USE_CUDA)
#  include <cusparse.h>
#elif defined(AMREX_USE_HIP)
#  include <rocsparse/rocsparse.h>
#endif

/**
 * \file AMReX_GpuSparse.H
 *
 * The cuSPARSE or rocSPARSE handles that AMReX's sparse linear algebra
 * (SpMV, SpGEMM, transpose) uses, for your own library calls.  AMReX keeps
 * them until amrex::Finalize.
 */

namespace amrex::Gpu {

#if defined(AMREX_USE_CUDA)
/**
 * \brief The cuSPARSE handle of the calling host thread for the current GPU
 *        stream (Gpu::gpuStream()).
 *
 * Each call binds the handle to that stream and sets the host pointer mode;
 * call it again after switching streams.  AMReX's sparse operations (SpMV,
 * SpGEMM, transpose) call it too, so set your own stream or pointer mode again
 * after each of them.  If the handle was bound to another stream (by you, or
 * for another external stream override), the call first synchronizes the
 * device, so it may block and must not run while a stream is captured into a
 * graph.  The handle is owned by AMReX: do not destroy it, and restore any
 * other setting you change.
 *
 * Call it while AMReX is initialized (amrex::ExecOnFinalize callbacks may),
 * from a thread that uses AMReX's GPU (Gpu::Device::deviceId()).  No other
 * thread may make sparse calls while amrex::Finalize runs, or while a thread
 * sets or resets an external stream (Gpu::ExternalGpuStreamRegion), which
 * applies to all threads.
 */
cusparseHandle_t getCusparseHandle ();
#elif defined(AMREX_USE_HIP)
/**
 * \brief The rocSPARSE handle of the calling host thread for the current GPU
 *        stream (Gpu::gpuStream()).
 *
 * Each call binds the handle to that stream and sets the host pointer mode;
 * call it again after switching streams.  AMReX's sparse operations (SpMV,
 * SpGEMM, transpose) call it too, so set your own stream or pointer mode again
 * after each of them.  If the handle was bound to another stream (by you, or
 * for another external stream override), the call first synchronizes the
 * device, so it may block and must not run while a stream is captured into a
 * graph.  The handle is owned by AMReX: do not destroy it, and restore any
 * other setting you change.
 *
 * Call it while AMReX is initialized (amrex::ExecOnFinalize callbacks may),
 * from a thread that uses AMReX's GPU (Gpu::Device::deviceId()).  No other
 * thread may make sparse calls while amrex::Finalize runs, or while a thread
 * sets or resets an external stream (Gpu::ExternalGpuStreamRegion), which
 * applies to all threads.
 */
rocsparse_handle getRocsparseHandle ();
#endif

}

#endif
