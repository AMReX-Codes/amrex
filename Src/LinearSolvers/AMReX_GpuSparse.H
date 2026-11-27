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
 * Cached handles of the vendor sparse-matrix libraries.
 *
 * Creating a cuSPARSE or rocSPARSE handle allocates device resources and
 * costs of the order of a millisecond, far more than a sparse matrix-vector
 * product on a coarse multigrid level.  The handles are therefore created
 * once per host thread, on first use, and destroyed in amrex::Finalize.
 */

namespace amrex::Gpu {

#if defined(AMREX_USE_CUDA)
/**
 * \brief The cuSPARSE handle of the calling host thread, bound to the
 *        current GPU stream (Gpu::gpuStream()).
 *
 * The handle is owned by AMReX: do not destroy it.
 */
cusparseHandle_t getCusparseHandle ();
#elif defined(AMREX_USE_HIP)
/**
 * \brief The rocSPARSE handle of the calling host thread, bound to the
 *        current GPU stream (Gpu::gpuStream()).
 *
 * The handle is owned by AMReX: do not destroy it.
 */
rocsparse_handle getRocsparseHandle ();
#endif

}

#endif
