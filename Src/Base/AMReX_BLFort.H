#ifndef _BL_BLFORT_H_
#define _BL_BLFORT_H_ 

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

#define BL_FORT_FAB_ARG(A) amrex_real* A##_fab, ARLIM_P(A##_lo), ARLIM_P(A##_hi)
#define BL_FORT_IFAB_ARG(A) int* A##_fab, ARLIM_P(A##_lo), ARLIM_P(A##_hi)
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
// FORT_LAUNCH(mfi, box, function, ...)
// where everything after the function name (represented here by ...)
// are the actual function arguments, as defined in the Fortran subroutine.
// mfi is the MFIter, and box is the current Box for this MFIter iteration.

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

#define FORT_LAUNCH(mfi, box, function, ...)  \
    dim3 numBlocks, numThreads; \
    Device::c_threads_and_blocks(box.loVect(), box.hiVect(), numBlocks, numThreads); \
    cudaStream_t stream = Device::stream_from_index(mfi.tileIndex()); \
    function<<<numBlocks, numThreads, 0, stream>>>(__VA_ARGS__)
#define FORT_LAUNCHABLE(func) __global__ func

#else

#define FORT_LAUNCH(mfi, box, function, ...) function(__VA_ARGS__)
#define FORT_LAUNCHABLE(func) func

#endif

#endif
