// Included first by every LOOM device pass (amrex-loom-cxx).
// Adapts GCC's libstdc++ to Loom's C++ frontend. Each line works around a
// frontend gap; none changes the meaning of code that is lowered to the GPU.
#ifndef AMREX_LOOM_DEVICE_PRELUDE_H_
#define AMREX_LOOM_DEVICE_PRELUDE_H_

// x86 spin-wait hint used by libstdc++ atomics on the host path.
#define __builtin_ia32_pause() ((void)0)
// <chrono> pulls in <format> via chrono_io; its consteval encoding checks are
// not evaluable by the frontend and AMReX does not format chrono types.
#define _GLIBCXX_CHRONO_IO_H 1
#include <bits/c++config.h>
// No __float128 type in the frontend.
#undef _GLIBCXX_USE_FLOAT128
// `if constexpr (trait{})` (class-type conversion) is not a frontend constant.
#undef _GLIBCXX17_CONSTEXPR
#define _GLIBCXX17_CONSTEXPR

// libstdc++ numeric_limits uses these builtins; spell them as constant bit
// patterns the importer lowers (IEEE binary32/binary64).
#define __builtin_huge_val() __builtin_bit_cast(double, 0x7ff0000000000000ull)
#define __builtin_huge_valf() __builtin_bit_cast(float, 0x7f800000u)
#define __builtin_inf() __builtin_bit_cast(double, 0x7ff0000000000000ull)
#define __builtin_inff() __builtin_bit_cast(float, 0x7f800000u)
#define __builtin_nan(s) __builtin_bit_cast(double, 0x7ff8000000000000ull)
#define __builtin_nanf(s) __builtin_bit_cast(float, 0x7fc00000u)
#define __builtin_nans(s) __builtin_bit_cast(double, 0x7ff4000000000000ull)
#define __builtin_nansf(s) __builtin_bit_cast(float, 0x7fa00000u)


// libm on the device: glibc declares these extern "C" (libstdc++ <cmath>
// re-exports them as std::). Define them first as inline wrappers over Loom
// scalar ops; the later glibc declarations are compatible redeclarations.
// The __builtin_ spellings <cmath> uses for float overloads map here too.
#include <loomcxx/scalar.h>
#define AMREX_LOOM_LIBM1(name, op)                                              \
  extern "C" inline double name(double x) noexcept { return ::loom::scalar::op(x); } \
  extern "C" inline float name##f(float x) noexcept { return ::loom::scalar::op(x); }
#define AMREX_LOOM_LIBM2(name, op)                                              \
  extern "C" inline double name(double x, double y) noexcept {                  \
    return ::loom::scalar::op(x, y);                                            \
  }                                                                             \
  extern "C" inline float name##f(float x, float y) noexcept {                  \
    return ::loom::scalar::op(x, y);                                            \
  }
AMREX_LOOM_LIBM1(sqrt, sqrtf)
AMREX_LOOM_LIBM1(exp, expf)
AMREX_LOOM_LIBM1(exp2, exp2f)
AMREX_LOOM_LIBM1(expm1, expm1f)
AMREX_LOOM_LIBM1(log, logf)
AMREX_LOOM_LIBM1(log2, log2f)
AMREX_LOOM_LIBM1(log10, log10f)
AMREX_LOOM_LIBM1(log1p, log1pf)
AMREX_LOOM_LIBM1(sin, sinf)
AMREX_LOOM_LIBM1(cos, cosf)
AMREX_LOOM_LIBM1(tan, tanf)
AMREX_LOOM_LIBM1(sinh, sinhf)
AMREX_LOOM_LIBM1(cosh, coshf)
AMREX_LOOM_LIBM1(tanh, tanhf)
AMREX_LOOM_LIBM1(asin, asinf)
AMREX_LOOM_LIBM1(acos, acosf)
AMREX_LOOM_LIBM1(atan, atanf)
AMREX_LOOM_LIBM1(asinh, asinhf)
AMREX_LOOM_LIBM1(acosh, acoshf)
AMREX_LOOM_LIBM1(atanh, atanhf)
AMREX_LOOM_LIBM1(cbrt, cbrtf)
AMREX_LOOM_LIBM1(floor, floorf)
AMREX_LOOM_LIBM1(ceil, ceilf)
AMREX_LOOM_LIBM1(trunc, truncf)
AMREX_LOOM_LIBM1(round, roundf)
AMREX_LOOM_LIBM1(erf, erff)
AMREX_LOOM_LIBM1(erfc, erfcf)
AMREX_LOOM_LIBM1(fabs, absf)
AMREX_LOOM_LIBM2(pow, powf)
AMREX_LOOM_LIBM2(atan2, atan2f)
AMREX_LOOM_LIBM2(copysign, copysignf)
AMREX_LOOM_LIBM2(fmax, maxnumf)
AMREX_LOOM_LIBM2(fmin, minnumf)
extern "C" inline double fma(double a, double b, double c) noexcept {
  return ::loom::scalar::fmaf(a, b, c);
}
extern "C" inline float fmaf(float a, float b, float c) noexcept {
  return ::loom::scalar::fmaf(a, b, c);
}
extern "C" inline int abs(int x) noexcept { return x < 0 ? -x : x; }
extern "C" inline long labs(long x) noexcept { return x < 0 ? -x : x; }
extern "C" inline long long llabs(long long x) noexcept { return x < 0 ? -x : x; }
#undef AMREX_LOOM_LIBM1
#undef AMREX_LOOM_LIBM2
#define __builtin_sqrtf ::sqrtf
#define __builtin_sqrt ::sqrt
#define __builtin_expf ::expf
#define __builtin_exp ::exp
#define __builtin_exp2f ::exp2f
#define __builtin_exp2 ::exp2
#define __builtin_expm1f ::expm1f
#define __builtin_expm1 ::expm1
#define __builtin_logf ::logf
#define __builtin_log ::log
#define __builtin_log2f ::log2f
#define __builtin_log2 ::log2
#define __builtin_log10f ::log10f
#define __builtin_log10 ::log10
#define __builtin_log1pf ::log1pf
#define __builtin_log1p ::log1p
#define __builtin_sinf ::sinf
#define __builtin_sin ::sin
#define __builtin_cosf ::cosf
#define __builtin_cos ::cos
#define __builtin_tanf ::tanf
#define __builtin_tan ::tan
#define __builtin_sinhf ::sinhf
#define __builtin_sinh ::sinh
#define __builtin_coshf ::coshf
#define __builtin_cosh ::cosh
#define __builtin_tanhf ::tanhf
#define __builtin_tanh ::tanh
#define __builtin_asinf ::asinf
#define __builtin_asin ::asin
#define __builtin_acosf ::acosf
#define __builtin_acos ::acos
#define __builtin_atanf ::atanf
#define __builtin_atan ::atan
#define __builtin_asinhf ::asinhf
#define __builtin_asinh ::asinh
#define __builtin_acoshf ::acoshf
#define __builtin_acosh ::acosh
#define __builtin_atanhf ::atanhf
#define __builtin_atanh ::atanh
#define __builtin_cbrtf ::cbrtf
#define __builtin_cbrt ::cbrt
#define __builtin_floorf ::floorf
#define __builtin_floor ::floor
#define __builtin_ceilf ::ceilf
#define __builtin_ceil ::ceil
#define __builtin_truncf ::truncf
#define __builtin_trunc ::trunc
#define __builtin_roundf ::roundf
#define __builtin_round ::round
#define __builtin_erff ::erff
#define __builtin_erf ::erf
#define __builtin_erfcf ::erfcf
#define __builtin_erfc ::erfc
#define __builtin_fabsf ::fabsf
#define __builtin_fabs ::fabs
#define __builtin_powf ::powf
#define __builtin_pow ::pow
#define __builtin_atan2f ::atan2f
#define __builtin_atan2 ::atan2
#define __builtin_copysignf ::copysignf
#define __builtin_copysign ::copysign
#define __builtin_fmaxf ::fmaxf
#define __builtin_fmax ::fmax
#define __builtin_fminf ::fminf
#define __builtin_fmin ::fmin
#define __builtin_fmaf ::fmaf
#define __builtin_fma ::fma


#endif
