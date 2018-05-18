#ifndef AMREX_BLFORT_H_
#define AMREX_BLFORT_H_

#include <AMReX_ArrayLim.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

#ifdef AMREX_USE_CUDA
#include <AMReX_Device.H>
#endif

#ifndef BL_LANG_FORT

#if   defined(BL_FORT_USE_UPPERCASE)
#define __BL_FORT_NAME__(A,B) A
#elif defined(BL_FORT_USE_LOWERCASE)
#define __BL_FORT_NAME__(A,B) B
#elif defined(BL_FORT_USE_UNDERSCORE)
#define __BL_FORT_NAME__(A,B) B##_
#else
#error "One of BL_FORT_USE_{UPPERCASE,LOWERCASE,UNDERSCORE} must be defined"
#endif

#ifdef __cplusplus
#define BL_FORT_PROC_DECL(A,B) extern "C" void __BL_FORT_NAME__(A,B)
#else
#define BL_FORT_PROC_DECL(A,B) void __BL_FORT_NAME__(A,B)
#endif
#define BL_FORT_PROC_CALL(A,B) __BL_FORT_NAME__(A,B)
#define BL_FORT_PROC_NAME(A,B) __BL_FORT_NAME__(A,B)

#define BL_FORT_FAB_ARG(A) amrex_real* A##_fab, AMREX_ARLIM_P(A##_lo), AMREX_ARLIM_P(A##_hi)
#define BL_FORT_IFAB_ARG(A) int* A##_fab, AMREX_ARLIM_P(A##_lo), AMREX_ARLIM_P(A##_hi)
#define BL_FORT_FAB_VAL_3D(A) A##_fab, A##_lo, A##_hi
#define BL_FORT_FAB_ARG_3D(A) amrex_real* A##_fab, const int* A##_lo, const int* A##_hi
#define BL_FORT_FAB_ARG_ANYD(A) amrex_real* A##_fab, const int* A##_lo, const int* A##_hi
#define BL_FORT_IFAB_ARG_3D(A) int* A##_fab, const int* A##_lo, const int* A##_hi
#define BL_FORT_IFAB_ARG_ANYD(A) int* A##_fab, const int* A##_lo, const int* A##_hi

#endif

// Define a macro that launches a Fortran function. For normal
// CPU code, it does nothing -- only launches the function. For
// GPU code, it launches a kernel using the first two arguments
// as bounds of the box to launch the kernel on.

// The form of the function launch MUST be:
// FORT_LAUNCH(box, function, ...)
// where everything after the function name (represented here by ...)
// are the actual function arguments, as defined in the Fortran subroutine.
// box is the current Box for this MFIter iteration.

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

#define FORT_LAUNCH(box, function, ...)  \
    { \
      dim3 numBlocks, numThreads; \
      Device::c_threads_and_blocks(box.loVect(), box.hiVect(), numBlocks, numThreads); \
      function<<<numBlocks, numThreads, 0, Device::cudaStream()>>>(__VA_ARGS__); \
      CudaErrorCheck(); \
    }

#define FORT_LAUNCHABLE(func) __global__ func

#else

#define FORT_LAUNCH(box, function, ...) function(__VA_ARGS__)
#define FORT_LAUNCHABLE(func) func

#endif

#endif

// This macro "deletes" a function when we're using CUDA,
// so that a script can replace it.

#ifdef AMREX_USE_CUDA
#define DEVICE_LAUNCHABLE(func) null##func
#define DEVICE_LAUNCHABLE_STRIDED(func) null##func
#else
#define DEVICE_LAUNCHABLE(func) func
#define DEVICE_LAUNCHABLE_STRIDED(func) func
#endif

// Decorate the calls to the GPU.

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

#if __CUDACC_VER_MAJOR__ >= 9 && __CUDACC_VER_MINOR__ >= 1

#define DEVICE_LAUNCH(function)                         \
      dim3 function_##numBlocks, function_##numThreads; \
      Device::grid_stride_threads_and_blocks(function_##numBlocks, function_##numThreads); \
      CudaAPICheck(cudaFuncSetAttribute(&cuda_##function, cudaFuncAttributePreferredSharedMemoryCarveout, 0)); \
      cuda_##function<<<function_##numBlocks, function_##numThreads, 0, Device::cudaStream()>>> \

#else

#define DEVICE_LAUNCH(function)                         \
      dim3 function_##numBlocks, function_##numThreads; \
      Device::grid_stride_threads_and_blocks(function_##numBlocks, function_##numThreads); \
      cuda_##function<<<function_##numBlocks, function_##numThreads, 0, Device::cudaStream()>>> \

#endif

#else

#define DEVICE_LAUNCH(function) function

#endif


