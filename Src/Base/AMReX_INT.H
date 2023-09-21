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

#endif
