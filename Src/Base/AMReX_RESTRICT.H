#ifndef AMREX_RESTRICT_H_
#define AMREX_RESTRICT_H_

#ifdef __cplusplus

#ifdef _WIN32
#define AMREX_RESTRICT __restrict
#else
#define AMREX_RESTRICT __restrict__
#endif

#else

#define AMREX_RESTRICT restrict

#endif


#if defined(__CUDA_ARCH__)
#define AMREX_PRAGMA_SIMD 

#elif defined(_OPENMP)
#define AMREX_PRAGMA_SIMD _Pragma("omp simd")

#elif defined(__INTEL_COMPILER)
#define AMREX_PRAGMA_SIMD _Pragma("simd")

#elif defined(_CRAYC)
#define AMREX_PRAGMA_SIMD _Pragma("ivdep")

#elif defined(__PGI)
#define AMREX_PRAGMA_SIMD _Pragma("vector")

#elif defined(__NEC__)
#define AMREX_PRAGMA_SIMD 

#elif defined(__ibmxl__)
#define AMREX_PRAGMA_SIMD _Pragma("simd_level(10)")

#elif defined(__llvm__)
#define AMREX_PRAGMA_SIMD _Pragma("clang loop vectorize(enable)")

#elif defined(__GNUC__)
#define AMREX_PRAGMA_SIMD _Pragma("GCC ivdep")

#else
#define AMREX_PRAGMA_SIMD 

#endif


#endif
