#ifndef AMREX_INT_H_
#define AMREX_INT_H_

#ifdef __cplusplus
#include <type_traits>
#include <cinttypes>
#else
#ifdef AMREX_TYPECHECK
#define __attribute__(x) 
#define __extension__ 
#define __restrict 
#endif
#include <inttypes.h>
#endif

#ifdef _WIN32
typedef          long long amrex_long;
typedef unsigned long long amrex_ulong;
#else
typedef          long amrex_long;
typedef unsigned long amrex_ulong;
#endif

#ifdef __cplusplus
namespace amrex {
    using Long = amrex_long;
    using ULong = amrex_ulong;
}
#endif

#endif
