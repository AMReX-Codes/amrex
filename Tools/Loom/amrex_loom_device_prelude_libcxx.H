// Device prelude for the LOOM backend with libc++ (single-source builds, where
// cxx compiles the host side too). amrex-loom-cxx includes it first in the
// device pass only (AMREX_LOOM_DEVICE_PRELUDE). It supplies what Loom does not
// yet provide for hosted C++: C math functions over loom::scalar, and the
// __builtin_* spellings libc++ uses. Remove it once Loom ships a <cmath>.
#ifndef AMREX_LOOM_DEVICE_PRELUDE_LIBCXX_H_
#define AMREX_LOOM_DEVICE_PRELUDE_LIBCXX_H_

// x86 spin-wait hint, a no-op on the device.
#define __builtin_ia32_pause() ((void)0)

// numeric_limits and <cmath> use these builtins; spell them as constant bit
// patterns the importer lowers (IEEE binary32/binary64). Removing them drops
// four kernels in AMReX_Interp_C.H.
#define __builtin_huge_val() __builtin_bit_cast(double, 0x7ff0000000000000ull)
#define __builtin_huge_valf() __builtin_bit_cast(float, 0x7f800000u)
#define __builtin_inf() __builtin_bit_cast(double, 0x7ff0000000000000ull)
#define __builtin_inff() __builtin_bit_cast(float, 0x7f800000u)
#define __builtin_nan(s) __builtin_bit_cast(double, 0x7ff8000000000000ull)
#define __builtin_nanf(s) __builtin_bit_cast(float, 0x7fc00000u)
#define __builtin_nans(s) __builtin_bit_cast(double, 0x7ff4000000000000ull)
#define __builtin_nansf(s) __builtin_bit_cast(float, 0x7fa00000u)


// libm on the device: glibc declares these extern "C" and libc++ <cmath>
// re-exports them as std::. Define them first as inline wrappers over Loom
// scalar ops; the later glibc declarations are compatible redeclarations.
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
// libc++ calls these builtins from inside namespace std, qualified or not.
// Each maps to a distinct C wrapper over the prelude's ::name. Mapping them to
// a plain name instead resolves to std:: overloads inside std (for example,
// double std::sqrt(double) would call itself or the float overload).
extern "C" inline float __amrex_loom_builtin_sqrtf(float a0) noexcept { return ::sqrtf(a0); }
extern "C" inline double __amrex_loom_builtin_sqrt(double a0) noexcept { return ::sqrt(a0); }
extern "C" inline float __amrex_loom_builtin_expf(float a0) noexcept { return ::expf(a0); }
extern "C" inline double __amrex_loom_builtin_exp(double a0) noexcept { return ::exp(a0); }
extern "C" inline float __amrex_loom_builtin_exp2f(float a0) noexcept { return ::exp2f(a0); }
extern "C" inline double __amrex_loom_builtin_exp2(double a0) noexcept { return ::exp2(a0); }
extern "C" inline float __amrex_loom_builtin_expm1f(float a0) noexcept { return ::expm1f(a0); }
extern "C" inline double __amrex_loom_builtin_expm1(double a0) noexcept { return ::expm1(a0); }
extern "C" inline float __amrex_loom_builtin_logf(float a0) noexcept { return ::logf(a0); }
extern "C" inline double __amrex_loom_builtin_log(double a0) noexcept { return ::log(a0); }
extern "C" inline float __amrex_loom_builtin_log2f(float a0) noexcept { return ::log2f(a0); }
extern "C" inline double __amrex_loom_builtin_log2(double a0) noexcept { return ::log2(a0); }
extern "C" inline float __amrex_loom_builtin_log10f(float a0) noexcept { return ::log10f(a0); }
extern "C" inline double __amrex_loom_builtin_log10(double a0) noexcept { return ::log10(a0); }
extern "C" inline float __amrex_loom_builtin_log1pf(float a0) noexcept { return ::log1pf(a0); }
extern "C" inline double __amrex_loom_builtin_log1p(double a0) noexcept { return ::log1p(a0); }
extern "C" inline float __amrex_loom_builtin_sinf(float a0) noexcept { return ::sinf(a0); }
extern "C" inline double __amrex_loom_builtin_sin(double a0) noexcept { return ::sin(a0); }
extern "C" inline float __amrex_loom_builtin_cosf(float a0) noexcept { return ::cosf(a0); }
extern "C" inline double __amrex_loom_builtin_cos(double a0) noexcept { return ::cos(a0); }
extern "C" inline float __amrex_loom_builtin_tanf(float a0) noexcept { return ::tanf(a0); }
extern "C" inline double __amrex_loom_builtin_tan(double a0) noexcept { return ::tan(a0); }
extern "C" inline float __amrex_loom_builtin_sinhf(float a0) noexcept { return ::sinhf(a0); }
extern "C" inline double __amrex_loom_builtin_sinh(double a0) noexcept { return ::sinh(a0); }
extern "C" inline float __amrex_loom_builtin_coshf(float a0) noexcept { return ::coshf(a0); }
extern "C" inline double __amrex_loom_builtin_cosh(double a0) noexcept { return ::cosh(a0); }
extern "C" inline float __amrex_loom_builtin_tanhf(float a0) noexcept { return ::tanhf(a0); }
extern "C" inline double __amrex_loom_builtin_tanh(double a0) noexcept { return ::tanh(a0); }
extern "C" inline float __amrex_loom_builtin_asinf(float a0) noexcept { return ::asinf(a0); }
extern "C" inline double __amrex_loom_builtin_asin(double a0) noexcept { return ::asin(a0); }
extern "C" inline float __amrex_loom_builtin_acosf(float a0) noexcept { return ::acosf(a0); }
extern "C" inline double __amrex_loom_builtin_acos(double a0) noexcept { return ::acos(a0); }
extern "C" inline float __amrex_loom_builtin_atanf(float a0) noexcept { return ::atanf(a0); }
extern "C" inline double __amrex_loom_builtin_atan(double a0) noexcept { return ::atan(a0); }
extern "C" inline float __amrex_loom_builtin_asinhf(float a0) noexcept { return ::asinhf(a0); }
extern "C" inline double __amrex_loom_builtin_asinh(double a0) noexcept { return ::asinh(a0); }
extern "C" inline float __amrex_loom_builtin_acoshf(float a0) noexcept { return ::acoshf(a0); }
extern "C" inline double __amrex_loom_builtin_acosh(double a0) noexcept { return ::acosh(a0); }
extern "C" inline float __amrex_loom_builtin_atanhf(float a0) noexcept { return ::atanhf(a0); }
extern "C" inline double __amrex_loom_builtin_atanh(double a0) noexcept { return ::atanh(a0); }
extern "C" inline float __amrex_loom_builtin_cbrtf(float a0) noexcept { return ::cbrtf(a0); }
extern "C" inline double __amrex_loom_builtin_cbrt(double a0) noexcept { return ::cbrt(a0); }
extern "C" inline float __amrex_loom_builtin_floorf(float a0) noexcept { return ::floorf(a0); }
extern "C" inline double __amrex_loom_builtin_floor(double a0) noexcept { return ::floor(a0); }
extern "C" inline float __amrex_loom_builtin_ceilf(float a0) noexcept { return ::ceilf(a0); }
extern "C" inline double __amrex_loom_builtin_ceil(double a0) noexcept { return ::ceil(a0); }
extern "C" inline float __amrex_loom_builtin_truncf(float a0) noexcept { return ::truncf(a0); }
extern "C" inline double __amrex_loom_builtin_trunc(double a0) noexcept { return ::trunc(a0); }
extern "C" inline float __amrex_loom_builtin_roundf(float a0) noexcept { return ::roundf(a0); }
extern "C" inline double __amrex_loom_builtin_round(double a0) noexcept { return ::round(a0); }
extern "C" inline float __amrex_loom_builtin_erff(float a0) noexcept { return ::erff(a0); }
extern "C" inline double __amrex_loom_builtin_erf(double a0) noexcept { return ::erf(a0); }
extern "C" inline float __amrex_loom_builtin_erfcf(float a0) noexcept { return ::erfcf(a0); }
extern "C" inline double __amrex_loom_builtin_erfc(double a0) noexcept { return ::erfc(a0); }
extern "C" inline float __amrex_loom_builtin_fabsf(float a0) noexcept { return ::fabsf(a0); }
extern "C" inline double __amrex_loom_builtin_fabs(double a0) noexcept { return ::fabs(a0); }
extern "C" inline float __amrex_loom_builtin_powf(float a0, float a1) noexcept { return ::powf(a0, a1); }
extern "C" inline double __amrex_loom_builtin_pow(double a0, double a1) noexcept { return ::pow(a0, a1); }
extern "C" inline float __amrex_loom_builtin_atan2f(float a0, float a1) noexcept { return ::atan2f(a0, a1); }
extern "C" inline double __amrex_loom_builtin_atan2(double a0, double a1) noexcept { return ::atan2(a0, a1); }
extern "C" inline float __amrex_loom_builtin_copysignf(float a0, float a1) noexcept { return ::copysignf(a0, a1); }
extern "C" inline double __amrex_loom_builtin_copysign(double a0, double a1) noexcept { return ::copysign(a0, a1); }
extern "C" inline float __amrex_loom_builtin_fmaxf(float a0, float a1) noexcept { return ::fmaxf(a0, a1); }
extern "C" inline double __amrex_loom_builtin_fmax(double a0, double a1) noexcept { return ::fmax(a0, a1); }
extern "C" inline float __amrex_loom_builtin_fminf(float a0, float a1) noexcept { return ::fminf(a0, a1); }
extern "C" inline double __amrex_loom_builtin_fmin(double a0, double a1) noexcept { return ::fmin(a0, a1); }
extern "C" inline float __amrex_loom_builtin_fmaf(float a0, float a1, float a2) noexcept { return ::fmaf(a0, a1, a2); }
extern "C" inline double __amrex_loom_builtin_fma(double a0, double a1, double a2) noexcept { return ::fma(a0, a1, a2); }
#define __builtin_sqrtf __amrex_loom_builtin_sqrtf
#define __builtin_sqrt __amrex_loom_builtin_sqrt
#define __builtin_expf __amrex_loom_builtin_expf
#define __builtin_exp __amrex_loom_builtin_exp
#define __builtin_exp2f __amrex_loom_builtin_exp2f
#define __builtin_exp2 __amrex_loom_builtin_exp2
#define __builtin_expm1f __amrex_loom_builtin_expm1f
#define __builtin_expm1 __amrex_loom_builtin_expm1
#define __builtin_logf __amrex_loom_builtin_logf
#define __builtin_log __amrex_loom_builtin_log
#define __builtin_log2f __amrex_loom_builtin_log2f
#define __builtin_log2 __amrex_loom_builtin_log2
#define __builtin_log10f __amrex_loom_builtin_log10f
#define __builtin_log10 __amrex_loom_builtin_log10
#define __builtin_log1pf __amrex_loom_builtin_log1pf
#define __builtin_log1p __amrex_loom_builtin_log1p
#define __builtin_sinf __amrex_loom_builtin_sinf
#define __builtin_sin __amrex_loom_builtin_sin
#define __builtin_cosf __amrex_loom_builtin_cosf
#define __builtin_cos __amrex_loom_builtin_cos
#define __builtin_tanf __amrex_loom_builtin_tanf
#define __builtin_tan __amrex_loom_builtin_tan
#define __builtin_sinhf __amrex_loom_builtin_sinhf
#define __builtin_sinh __amrex_loom_builtin_sinh
#define __builtin_coshf __amrex_loom_builtin_coshf
#define __builtin_cosh __amrex_loom_builtin_cosh
#define __builtin_tanhf __amrex_loom_builtin_tanhf
#define __builtin_tanh __amrex_loom_builtin_tanh
#define __builtin_asinf __amrex_loom_builtin_asinf
#define __builtin_asin __amrex_loom_builtin_asin
#define __builtin_acosf __amrex_loom_builtin_acosf
#define __builtin_acos __amrex_loom_builtin_acos
#define __builtin_atanf __amrex_loom_builtin_atanf
#define __builtin_atan __amrex_loom_builtin_atan
#define __builtin_asinhf __amrex_loom_builtin_asinhf
#define __builtin_asinh __amrex_loom_builtin_asinh
#define __builtin_acoshf __amrex_loom_builtin_acoshf
#define __builtin_acosh __amrex_loom_builtin_acosh
#define __builtin_atanhf __amrex_loom_builtin_atanhf
#define __builtin_atanh __amrex_loom_builtin_atanh
#define __builtin_cbrtf __amrex_loom_builtin_cbrtf
#define __builtin_cbrt __amrex_loom_builtin_cbrt
#define __builtin_floorf __amrex_loom_builtin_floorf
#define __builtin_floor __amrex_loom_builtin_floor
#define __builtin_ceilf __amrex_loom_builtin_ceilf
#define __builtin_ceil __amrex_loom_builtin_ceil
#define __builtin_truncf __amrex_loom_builtin_truncf
#define __builtin_trunc __amrex_loom_builtin_trunc
#define __builtin_roundf __amrex_loom_builtin_roundf
#define __builtin_round __amrex_loom_builtin_round
#define __builtin_erff __amrex_loom_builtin_erff
#define __builtin_erf __amrex_loom_builtin_erf
#define __builtin_erfcf __amrex_loom_builtin_erfcf
#define __builtin_erfc __amrex_loom_builtin_erfc
#define __builtin_fabsf __amrex_loom_builtin_fabsf
#define __builtin_fabs __amrex_loom_builtin_fabs
#define __builtin_powf __amrex_loom_builtin_powf
#define __builtin_pow __amrex_loom_builtin_pow
#define __builtin_atan2f __amrex_loom_builtin_atan2f
#define __builtin_atan2 __amrex_loom_builtin_atan2
#define __builtin_copysignf __amrex_loom_builtin_copysignf
#define __builtin_copysign __amrex_loom_builtin_copysign
#define __builtin_fmaxf __amrex_loom_builtin_fmaxf
#define __builtin_fmax __amrex_loom_builtin_fmax
#define __builtin_fminf __amrex_loom_builtin_fminf
#define __builtin_fmin __amrex_loom_builtin_fmin
#define __builtin_fmaf __amrex_loom_builtin_fmaf
#define __builtin_fma __amrex_loom_builtin_fma

#endif
