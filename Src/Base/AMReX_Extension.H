#ifndef AMREX_EXTENSION_H_
#define AMREX_EXTENSION_H_
#include <AMReX_Config.H>

#if !defined(BL_LANG_FORT)

// restrict

#ifdef __cplusplus


#if defined(__INTEL_COMPILER)
#define AMREX_CXX_INTEL
#elif defined(_CRAYC)
#define AMREX_CXX_CRAY
#elif defined(__PGI)
#define AMREX_CXX_PGI
#elif defined(__NVCOMPILER)
#define AMREX_CXX_NVHPC
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

#elif defined(__CUDA_ARCH__) && defined(AMREX_USE_CUDA)
#define AMREX_PRAGMA_SIMD

#elif defined(__HIP_DEVICE_COMPILE__)
#define AMREX_PRAGMA_SIMD

//#elif defined(AMREX_USE_OMP) && defined(_OPENMP) && (_OPENMP >= 201307) && !defined(__PGI)
//#define AMREX_PRAGMA_SIMD _Pragma("omp simd")

#elif defined(__INTEL_COMPILER)
#define AMREX_PRAGMA_SIMD _Pragma("ivdep")

#elif defined(_CRAYC) || defined(__cray__)
#define AMREX_PRAGMA_SIMD _Pragma("ivdep")

#elif defined(__PGI)
#define AMREX_PRAGMA_SIMD _Pragma("loop ivdep")

#elif defined(__NVCOMPILER)
#define AMREX_PRAGMA_SIMD _Pragma("loop ivdep")

#elif defined(__NEC__)
#define AMREX_PRAGMA_SIMD

#elif defined(__ibmxl__)
#define AMREX_PRAGMA_SIMD _Pragma("ibm independent_loop")

#elif defined(__clang__)
#define AMREX_PRAGMA_SIMD

#elif defined(__GNUC__)
#define AMREX_PRAGMA_SIMD _Pragma("GCC ivdep")

#else
#define AMREX_PRAGMA_SIMD

#endif /* simd */

// force inline
#if defined(__CUDA_ARCH__) && defined(AMREX_USE_CUDA)
#define AMREX_FORCE_INLINE __forceinline__

#elif defined(__HIP_DEVICE_COMPILE__)
#define AMREX_FORCE_INLINE __forceinline__

#elif defined(__INTEL_COMPILER)
#define AMREX_FORCE_INLINE inline __attribute__((always_inline))

#elif defined(_CRAYC)
#define AMREX_FORCE_INLINE inline

#elif defined(__PGI)
#define AMREX_FORCE_INLINE inline

#elif defined(__NVCOMPILER)
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
#if defined(_MSC_VER)
#define AMREX_NO_INLINE __declspec(noinline)
#elif (defined(__GNUC__) || defined(__clang__) || defined(__CUDACC__) || defined(__HIP__) || defined(__INTEL_CLANG_COMPILER))
#define AMREX_NO_INLINE __attribute__((noinline))
#else
#define AMREX_NO_INLINE
#endif

// __attribute__((weak))

#if defined(AMREX_TYPECHECK)
#define AMREX_ATTRIBUTE_WEAK
#elif defined(_WIN32)
#define AMREX_ATTRIBUTE_WEAK
#elif defined(__clang__) && defined(__apple_build_version__)
#define AMREX_ATTRIBUTE_WEAK __attribute__((weak_import))
#else
#define AMREX_ATTRIBUTE_WEAK __attribute__((weak))
#endif

// public globals
// https://stackoverflow.com/questions/54560832/cmake-windows-export-all-symbols-does-not-cover-global-variables/54568678#54568678
#if defined(_MSC_VER)
#   if defined(AMREX_IS_DLL)
#       if defined(AMREX_IS_DLL_BUILDING)
#           define AMREX_EXPORT __declspec(dllexport)
#       else
#           define AMREX_EXPORT __declspec(dllimport)
#       endif
#   else
#       define AMREX_EXPORT
#   endif
#else
#   define AMREX_EXPORT
#endif

#if defined(__cplusplus) && defined(_WIN32)
#include <ciso646>
#endif

#if defined(__INTEL_COMPILER) && defined(__EDG__) && (__cplusplus < 201703L)
// Classical EDG based Intel compiler does not support fallthrough when std=c++14
#    define AMREX_FALLTHROUGH ((void)0)
#elif defined(__has_cpp_attribute) && __has_cpp_attribute(fallthrough) >= 201603L
#    define AMREX_FALLTHROUGH [[fallthrough]]
#elif defined(__clang__)
#    define AMREX_FALLTHROUGH [[clang::fallthrough]]
#elif defined(__GNUC__) && (__GNUC__ >= 7)
#    define AMREX_FALLTHROUGH [[gnu::fallthrough]]
#else
#    define AMREX_FALLTHROUGH ((void)0)
#endif

// CI uses -Werror -Wc++17-extension, thus we need to add the __cplusplus clause
#if !defined(AMREX_NO_NODISCARD) && defined(__has_cpp_attribute) && __has_cpp_attribute(nodiscard) >= 201603L
#   define AMREX_NODISCARD [[nodiscard]]
#else
#   define AMREX_NODISCARD
#endif

// Note: following compilers support [[no_unique_address]]
//   - Clang >= 9.0
//   - GCC >= 9.0
//   - MSVC >= 19.26
// Using no unique address makes empty base class optimization for multiple policies much easier
#if !defined(AMREX_NO_NO_UNIQUE_ADDRESS) && defined(__has_cpp_attribute) && __has_cpp_attribute(no_unique_address) >= 201803L
#   define AMREX_NO_UNIQUE_ADDRESS [[no_unique_address]]
#   define AMREX_HAS_NO_UNIQUE_ADDRESS 1
#else
#   define AMREX_NO_UNIQUE_ADDRESS
#endif

#if defined(__cpp_if_constexpr) && __cpp_if_constexpr >= 201606L
#   define AMREX_IF_CONSTEXPR if constexpr
#else
#   define AMREX_IF_CONSTEXPR if
#endif

#if !defined(AMREX_NO_BUILTIN_CLZ)
#   if defined(__clang__) || defined(__GNUC__)
#      define AMREX_HAS_BUILTIN_CLZ 1
#   endif
#endif


#endif /* !BL_LANG_FORT */

#endif
