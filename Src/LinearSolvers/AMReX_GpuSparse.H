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
 * Handles of the vendor sparse-matrix libraries that AMReX's sparse linear
 * algebra uses: one per host thread and GPU stream, created on first use.
 * When a thread exits, its handles pass to the next threads that need one for
 * the same stream, and amrex::Finalize releases them all.
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
 * after each of them; if the handle is bound to another stream, the call
 * first waits for the device.  The handle is owned by AMReX: do not destroy
 * it, and restore any other setting you change.
 *
 * Call it while AMReX is initialized (amrex::ExecOnFinalize callbacks may),
 * from a thread that uses AMReX's GPU (Gpu::Device::deviceId()).
 * amrex::Finalize releases the handles of all threads, so no other thread may
 * make sparse calls while it runs.
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
 * after each of them; if the handle is bound to another stream, the call
 * first waits for the device.  The handle is owned by AMReX: do not destroy
 * it, and restore any other setting you change.
 *
 * Call it while AMReX is initialized (amrex::ExecOnFinalize callbacks may),
 * from a thread that uses AMReX's GPU (Gpu::Device::deviceId()).
 * amrex::Finalize releases the handles of all threads, so no other thread may
 * make sparse calls while it runs.
 */
rocsparse_handle getRocsparseHandle ();
#endif

}

#endif
