
#include <kdecl.H>

using namespace amrex;

#define SIMDORNOT(f) f##_simd

#include <kc.H>


#undef SIMDORNOT
#define SIMDORNOT(f) f##_nosimd
#undef AMREX_PRAGMA_SIMD
#define AMREX_PRAGMA_SIMD 

#include <kc.H>

