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
 * algebra uses: one per host thread, created on first use and released when
 * the thread exits or in amrex::Finalize (that of the main thread, which
 * loaded AMReX, in amrex::Finalize only).
 */

namespace amrex::Gpu {

#if defined(AMREX_USE_CUDA)
/**
 * \brief The cuSPARSE handle of the calling host thread.
 *
 * Each call binds the handle to the current GPU stream (Gpu::gpuStream()) and
 * sets the host pointer mode; call it again after switching streams.  AMReX's
 * sparse operations (SpMV, SpGEMM, transpose) call it too, so set your own
 * stream or pointer mode again after each of them.  The handle is owned by
 * AMReX: do not destroy it, and restore any other setting you change.
 *
 * The calling thread must use AMReX's GPU (Gpu::Device::deviceId()).
 * amrex::Finalize releases the handles of all threads, so no other thread may
 * make sparse calls while it runs.
 */
cusparseHandle_t getCusparseHandle ();
#elif defined(AMREX_USE_HIP)
/**
 * \brief The rocSPARSE handle of the calling host thread.
 *
 * Each call binds the handle to the current GPU stream (Gpu::gpuStream()) and
 * sets the host pointer mode; call it again after switching streams.  AMReX's
 * sparse operations (SpMV, SpGEMM, transpose) call it too, so set your own
 * stream or pointer mode again after each of them.  The handle is owned by
 * AMReX: do not destroy it, and restore any other setting you change.
 *
 * The calling thread must use AMReX's GPU (Gpu::Device::deviceId()).
 * amrex::Finalize releases the handles of all threads, so no other thread may
 * make sparse calls while it runs.
 */
rocsparse_handle getRocsparseHandle ();
#endif

}

#endif
