#ifndef AMREX_CONFIG_3D_H_
#define AMREX_CONFIG_3D_H_

/* #undef AMREX_XSDK */
/* #undef AMREX_DEBUG */
/* #undef AMREX_PROFILING */
/* #undef BL_PROFILING */
/* #undef AMREX_TRACE_PROFILING */
/* #undef BL_TRACE_PROFILING */
/* #undef AMREX_COMM_PROFILING */
/* #undef BL_COMM_PROFILING */
/* #undef AMREX_TINY_PROFILING */
/* #undef AMREX_USE_ROCTX */
/* #undef AMREX_MEM_PROFILING */
/* #undef AMREX_TESTING */
/* #undef AMREX_USE_MPI */
/* #undef BL_USE_MPI */
/* #undef AMREX_MPI_THREAD_MULTIPLE */
/* #undef AMREX_USE_SIMD */
/* #undef AMREX_USE_OMP */
/* #undef BL_USE_OMP */
/* #undef AMREX_USE_SYCL */
/* #undef AMREX_USE_DPCPP */
/* #undef AMREX_SYCL_SUB_GROUP_SIZE */
/* #undef AMREX_USE_ONEDPL */
/* #undef AMREX_USE_FLOAT */
/* #undef AMREX_SINGLE_PRECISION_PARTICLES */
/* #undef BL_USE_FLOAT */
#ifndef AMREX_SPACEDIM
#define AMREX_SPACEDIM 3
#endif
#ifndef BL_SPACEDIM
#define BL_SPACEDIM 3
#endif
/* #undef AMREX_USE_ASSERTION */
/* #undef AMREX_USE_FLATTEN_FOR */
/* #undef AMREX_BOUND_CHECK */
/* #undef AMREX_EXPORT_DYNAMIC */
/* #undef AMREX_FASTMATH */
/* #undef BL_FORT_USE_UNDERSCORE */
/* #undef BL_FORT_USE_LOWERCASE */
/* #undef BL_FORT_USE_UPPERCASE */
#define BL_NO_FORT
/* #undef AMREX_USE_FFT */
/* #undef AMREX_USE_MULTI_THREADED_FFTW */
/* #undef AMREX_USE_SENSEI_INSITU */
/* #undef AMREX_NO_SENSEI_AMR_INST */
/* #undef AMREX_USE_CONDUIT */
/* #undef AMREX_USE_CATALYST */
/* #undef AMREX_USE_ASCENT */
/* #undef AMREX_USE_EB */
/* #undef AMREX_USE_CUDA */
/* #undef AMREX_USE_HIP */
/* #undef AMREX_AMDGCN_WAVEFRONT_SIZE */
/* #undef AMREX_USE_NVML */
/* #undef AMREX_GPU_MAX_THREADS */
#ifndef AMREX_GPU_MAX_THREADS
#define AMREX_GPU_MAX_THREADS 0
#endif
/* #undef AMREX_USE_ACC */
/* #undef AMREX_USE_GPU */
/* #undef BL_COALESCE_FABS */
/* #undef AMREX_USE_GPU_RDC */
#define AMREX_PARTICLES
/* #undef AMREX_USE_HDF5 */
/* #undef AMREX_USE_HDF5_ASYNC */
/* #undef AMREX_USE_HDF5_ZFP */
/* #undef AMREX_USE_HYPRE */
/* #undef AMREX_USE_PETSC */
/* #undef AMREX_USE_SUNDIALS */
/* #undef AMREX_NO_PROBINIT */
/* #undef AMREX_IS_DLL */
#ifdef __cplusplus
#ifndef __GNUC__
static_assert(false,"libamrex was built with GNU. To avoid this error, reconfigure with -DAMReX_DIFFERENT_COMPILER=ON");
#endif

#endif
/*
  hipcc does not necessarily set the _OPENMP macro
  https://rocmdocs.amd.com/en/latest/Programming_Guides/HIP-FAQ.html?highlight=_openmp#openmp-is-undefined-when-compiling-with-fopenmp
 */
#if defined(AMREX_USE_OMP) && !defined(_OPENMP) && !defined(AMREX_USE_HIP)
#error "libamrex was built with OpenMP, so the downstream project must activate it, too"
#endif

#endif /* AMREX_CONFIG_3D_H_ */
