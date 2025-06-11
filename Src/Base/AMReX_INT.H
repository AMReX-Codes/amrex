#ifndef AMREX_INT_H_
#define AMREX_INT_H_
#include <AMReX_Config.H>

#ifdef __cplusplus
#include <type_traits>
#include <cinttypes>
#include <climits>
#else
#ifdef AMREX_TYPECHECK
#define __attribute__(x)
#define __extension__
#define __restrict
#endif
#include <inttypes.h>
#include <limits.h>
#endif

#if (INT_MAX == LONG_MAX)
typedef          long long amrex_long;  // NOLINT(modernize-use-using) // so that this can be a C header
typedef unsigned long long amrex_ulong; // NOLINT(modernize-use-using)
#else
typedef          long amrex_long;       // NOLINT(modernize-use-using)
typedef unsigned long amrex_ulong;      // NOLINT(modernize-use-using)
#endif

#ifdef __cplusplus
namespace amrex {
    using Long = amrex_long;
    using ULong = amrex_ulong;
}
#endif

#if (defined(__x86_64) || defined (__aarch64__)) && !defined(_WIN32) && (defined(__GNUC__) || defined(__clang__)) && !defined(__NVCOMPILER)

#define AMREX_INT128_SUPPORTED 1

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif

typedef unsigned __int128 amrex_uint128_t; // NOLINT(modernize-use-using)
typedef          __int128 amrex_int128_t;  // NOLINT(modernize-use-using)

#ifdef __cplusplus
namespace amrex {
    using UInt128_t = amrex_uint128_t;
    using  Int128_t = amrex_int128_t;
}
#endif

#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#endif /* (defined(__x86_64) || defined (__aarch64__)) && !defined(_WIN32) && (defined(__GNUC__) || defined(__clang__)) */

#endif
