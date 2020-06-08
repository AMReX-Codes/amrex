#ifndef AMREX_EXTENSION_H_
#define AMREX_EXTENSION_H_

#if !defined(BL_LANG_FORT)

// restrict

#ifdef __cplusplus


#if defined(__INTEL_COMPILER)
#define AMREX_CXX_INTEL
#elif defined(_CRAYC)
#define AMREX_CXX_CRAY
#elif defined(__PGI)
#define AMREX_CXX_PGI
#elif defined(__NEC__)
#define AMREX_CXX_NEC
#elif defined(__ibmxl__)
#define AMREX_CXX_IBM
#elif defined(__clang__)
#define AMREX_CXX_CLANG
#elif defined(__GNUC__)
#define AMREX_CXX_GCC
#endif


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

#elif defined(__HIP_DEVICE_COMPILE__)
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

#elif defined(__HIP_DEVICE_COMPILE__)
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

// no inline
#if defined(__GNUC__)
#define AMREX_NO_INLINE __attribute__((noinline))
#else
#define AMREX_NO_INLINE
#endif

// __attribute__((weak))

#if defined(AMREX_TYPECHECK)
#define AMREX_ATTRIBUTE_WEAK 
#elif defined(_WIN32)
#define AMREX_ATTRIBUTE_WEAK
#else
#define AMREX_ATTRIBUTE_WEAK __attribute__((weak))
#endif

#if defined(__cplusplus) && defined(_WIN32)
#include <iso646.h>
#endif

#if (__cplusplus >= 201703L)
#define AMREX_FALLTHROUGH [[fallthrough]]
#elif defined(__clang__)
#define AMREX_FALLTHROUGH [[clang::fallthrough]]
#elif defined(__GNUC__) && (__GNUC__ >= 7)
#define AMREX_FALLTHROUGH [[gnu::fallthrough]]
#else
#define AMREX_FALLTHROUGH ((void)0)
#endif


#endif /* !BL_LANG_FORT */

#endif
