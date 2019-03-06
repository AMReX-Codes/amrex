#ifndef AMREX_EXTENSION_H_
#define AMREX_EXTENSION_H_

#if !defined(BL_LANG_FORT)

// restrict

#ifdef __cplusplus

#ifdef _WIN32
#define AMREX_RESTRICT __restrict
#else
#define AMREX_RESTRICT __restrict__
#endif

#else

#define AMREX_RESTRICT restrict

#endif

// simd

#if defined(AMREX_DEBUG)
#define AMREX_PRAGMA_SIMD 

#elif defined(__CUDA_ARCH__)
#define AMREX_PRAGMA_SIMD 

//#elif defined(_OPENMP) && (_OPENMP >= 201307) && !defined(__PGI)
//#define AMREX_PRAGMA_SIMD _Pragma("omp simd")

#elif defined(__INTEL_COMPILER)
#define AMREX_PRAGMA_SIMD _Pragma("ivdep")

#elif defined(_CRAYC)
#define AMREX_PRAGMA_SIMD _Pragma("ivdep")

#elif defined(__PGI)
#define AMREX_PRAGMA_SIMD _Pragma("loop ivdep")

#elif defined(__NEC__)
#define AMREX_PRAGMA_SIMD 

#elif defined(__ibmxl__)
#define AMREX_PRAGMA_SIMD _Pragma("ibm independent_loop")

#elif defined(__clang__)
#define AMREX_PRAGMA_SIMD _Pragma("clang loop vectorize(enable)")

#elif defined(__GNUC__)
#define AMREX_PRAGMA_SIMD _Pragma("GCC ivdep")

#else
#define AMREX_PRAGMA_SIMD 

#endif /* simd */

// force inline
#if defined(__CUDA_ARCH__)
#define AMREX_FORCE_INLINE __forceinline__

#elif defined(__INTEL_COMPILER)
#define AMREX_FORCE_INLINE inline __attribute__((always_inline))

#elif defined(_CRAYC)
#define AMREX_FORCE_INLINE inline

#elif defined(__PGI)
#define AMREX_FORCE_INLINE inline

#elif defined(__NEC__)
#define AMREX_FORCE_INLINE inline

#elif defined(__ibmxl__)
#define AMREX_FORCE_INLINE inline __attribute__((always_inline))

#elif defined(__clang__)
#define AMREX_FORCE_INLINE inline __attribute__((always_inline))

#elif defined(__GNUC__)
#define AMREX_FORCE_INLINE inline __attribute__((always_inline))

#else
#define AMREX_FORCE_INLINE inline

#endif /* force inline */


#ifdef AMREX_USE_FORCE_INLINE
#define AMREX_INLINE AMREX_FORCE_INLINE
#else
#define AMREX_INLINE inline
#endif


#endif /* !BL_LANG_FORT */

#endif
