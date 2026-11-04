#ifndef AMREX_FPMP_H_
#define AMREX_FPMP_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_FPMP.H
 *
 * Support for the multi-precision floating-point types of NVIDIA CCCL
 * (<cuda/fpmp>), used when AMReX is built with AMReX_PRECISION=FP32MP2.
 * In that mode, amrex::Real and amrex::ParticleReal are cuda::fp32mp2, a
 * pair of floats (hi + lo) with a 46-bit significand and the exponent range
 * of float.
 *
 * This header is included by AMReX_REAL.H in every build, so the traits below
 * can be used in generic code independently of the precision AMReX was built
 * with.
 */

#if defined(AMREX_USE_FP32MP2) || defined(AMREX_FP32MP2_PARTICLES)
// Check the version first: an older CCCL (e.g., the copy bundled with the
// CUDA toolkit) found earlier on the include path would otherwise fail with
// obscure errors inside <cuda/fpmp>.
#  include <cuda/std/version>
#  if !defined(CCCL_VERSION) || (CCCL_VERSION < 3006000)
#    error "AMReX_PRECISION=FP32MP2 requires CCCL >= 3.6 (<cuda/fpmp>) ahead of any older CCCL on the include path"
#  endif
// AMReX uses CCCL's implicit conversion mode, so that code written for a
// built-in Real compiles unchanged, e.g., `Real x = 0.5;` or `mf.setVal(0.0);`.
// The setting changes the declarations in <cuda/fpmp> and so must be the same
// in every translation unit. CMake builds set it on the AMReX target (PUBLIC);
// this is the fallback for other build systems.
#  if !defined(CCCL_FPMP_EXPLICIT_CASTS)
#    define CCCL_FPMP_EXPLICIT_CASTS 0
#  elif CCCL_FPMP_EXPLICIT_CASTS != 0
#    error "AMReX_PRECISION=FP32MP2 requires CCCL_FPMP_EXPLICIT_CASTS=0: compile with -DCCCL_FPMP_EXPLICIT_CASTS=0, or include AMReX headers before <cuda/fpmp>"
#  endif
// The pair arithmetic relies on IEEE rounding of each float operation. Host
// fast-math lets the compiler reassociate it, which silently reduces the pairs
// to float accuracy. nvcc's --use_fast_math does not define these macros.
#  if defined(__FAST_MATH__) || defined(_M_FP_FAST)
#    error "AMReX_PRECISION=FP32MP2 cannot be compiled with host fast-math (-ffast-math, -Ofast, /fp:fast): it reduces cuda::fp32mp2 to float accuracy"
#  endif
#  include <cuda/fpmp>
#  include <iosfwd>
#  include <limits>
#  include <numbers>
#endif

#include <concepts>
#include <type_traits>

#if defined(AMREX_USE_FP32MP2) || defined(AMREX_FP32MP2_PARTICLES)
/**
 * CCCL only specializes cuda::std::numeric_limits for the fpmp2 types, while
 * AMReX and application codes query std::numeric_limits<amrex::Real>.
 * Specializing std::numeric_limits for a program-defined type is permitted by
 * the standard. The members are inherited and usable in device code.
 */
namespace std { // NOLINT(cert-dcl58-cpp)
template <cuda::fpmp2_accuracy A>
class numeric_limits<cuda::fpmp2<float, A>>
    : public cuda::std::numeric_limits<cuda::fpmp2<float, A>>
{};
} // namespace std

/**
 * The mathematical constants of <numbers> (e.g., std::numbers::pi_v<Real>),
 * which the standard permits specializing for a program-defined type.  The
 * values are rounded from double, which has more significand bits.  In
 * device code, copy them into a local variable (e.g.,
 * `constexpr Real pi = std::numbers::pi_v<Real>;`) or use
 * amrex::Math::pi<Real>(); using them directly in an expression with an
 * operator does not compile.
 */
namespace std::numbers { // NOLINT(cert-dcl58-cpp)
#define AMREX_FPMP2_NUMBERS_CONSTANT(name)                                  \
    template <cuda::fpmp2_accuracy A>                                      \
    inline constexpr cuda::fpmp2<float, A> name##_v<cuda::fpmp2<float, A>> \
        = cuda::fpmp2<float, A>(std::numbers::name);
AMREX_FPMP2_NUMBERS_CONSTANT(e)
AMREX_FPMP2_NUMBERS_CONSTANT(log2e)
AMREX_FPMP2_NUMBERS_CONSTANT(log10e)
AMREX_FPMP2_NUMBERS_CONSTANT(pi)
AMREX_FPMP2_NUMBERS_CONSTANT(inv_pi)
AMREX_FPMP2_NUMBERS_CONSTANT(inv_sqrtpi)
AMREX_FPMP2_NUMBERS_CONSTANT(ln2)
AMREX_FPMP2_NUMBERS_CONSTANT(ln10)
AMREX_FPMP2_NUMBERS_CONSTANT(sqrt2)
AMREX_FPMP2_NUMBERS_CONSTANT(sqrt3)
AMREX_FPMP2_NUMBERS_CONSTANT(inv_sqrt3)
AMREX_FPMP2_NUMBERS_CONSTANT(egamma)
AMREX_FPMP2_NUMBERS_CONSTANT(phi)
#undef AMREX_FPMP2_NUMBERS_CONSTANT
} // namespace std::numbers
#endif

namespace amrex {

/**
 * \brief True if T is one of the CCCL multi-precision pair types
 * (cuda::fpmp2<...>, e.g. cuda::fp32mp2).
 */
template <typename T>
struct IsFPMP2 : std::false_type {};

#if defined(AMREX_USE_FP32MP2) || defined(AMREX_FP32MP2_PARTICLES)
template <typename F, cuda::fpmp2_accuracy A>
struct IsFPMP2<cuda::fpmp2<F, A>> : std::true_type {};
#endif

template <typename T>
inline constexpr bool IsFPMP2_v = IsFPMP2<std::remove_cv_t<T>>::value;

/**
 * \brief A floating-point type that amrex::Real can be: float, double,
 * long double, or a CCCL multi-precision pair type (FP32MP2 builds).
 *
 * Use this instead of std::floating_point to constrain templates that
 * must accept amrex::Real in every precision mode.
 */
template <typename T>
concept FloatingPointLike = std::floating_point<T> || IsFPMP2_v<T>;

#if defined(AMREX_USE_FP32MP2) || defined(AMREX_FP32MP2_PARTICLES)
/**
 * \brief Read a cuda::fp32mp2 from a stream.
 *
 * The value is read as a double, which has more significand bits (53) than
 * cuda::fp32mp2 (46), and then split into the pair.  CCCL does not provide
 * stream operators.  Writing works without one, through the implicit
 * conversion of cuda::fp32mp2 to double.
 *
 * This operator is declared in namespace amrex, so it is found by ordinary
 * name lookup in AMReX code.  It is not found by argument-dependent lookup;
 * code outside namespace amrex needs `using amrex::operator>>;`.
 */
std::istream& operator>> (std::istream& is, cuda::fp32mp2& x);
#endif

} // namespace amrex

#endif
