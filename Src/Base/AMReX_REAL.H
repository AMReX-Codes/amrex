#ifndef AMREX_REAL_H_
#define AMREX_REAL_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>

#ifdef BL_USE_FLOAT
#  undef BL_USE_DOUBLE
#  undef BL_USE_FLOAT
/*
  The macro BL_USE_FLOAT indicates that C++ floating-point calculations should
  use "float" variables and Fortran floating-point calculations should use
  "real*4" variables.  One of BL_USE_FLOAT or BL_USE_DOUBLE must always be
  defined when compiling and using AMReX.
*/
#  define BL_USE_FLOAT 1
#else
#  undef  BL_USE_FLOAT
#  undef  BL_USE_DOUBLE
/*
  The macro BL_USE_DOUBLE indicates that C++ floating-point calculations
  should use "double" variables and Fortran floating-point calculations
  should use "real*8" variables.  One of BL_USE_FLOAT or BL_USE_DOUBLE must
  always be defined when compiling and using AMReX.
*/
#  define BL_USE_DOUBLE 1
#endif

#if !defined(BL_LANG_FORT)
#ifdef __cplusplus
#include <cfloat>
#else
#include <float.h>
#endif

/*
  Real is a typedef specifying the precision of the floating-point
  calculations in C++ code.  It will be either float or double
  depending upon which of the macros BL_USE_FLOAT or
  BL_USE_DOUBLE, respectively, is defined during compilations.  For
  portability, you should write floating-point code in terms of this
  typedef, instead of using float or double directly.

  Note that exactly one of these macros must be defined
  when compiling any module that uses floating-point.
*/

#ifdef BL_USE_FLOAT
typedef float amrex_real; // NOLINT(modernize-use-using)
// We need to define these to get around a CUDA 9.2 bug that breaks std::numeric_limits
#define AMREX_REAL_MIN     FLT_MIN
#define AMREX_REAL_MAX     FLT_MAX
#define AMREX_REAL_LOWEST (-FLT_MAX)
#else
typedef double amrex_real; // NOLINT(modernize-use-using)
#define AMREX_REAL_MIN     DBL_MIN
#define AMREX_REAL_MAX     DBL_MAX
#define AMREX_REAL_LOWEST (-DBL_MAX)
#endif

#ifdef AMREX_SINGLE_PRECISION_PARTICLES
typedef float amrex_particle_real; // NOLINT(modernize-use-using)
#else
typedef double amrex_particle_real; // NOLINT(modernize-use-using)
#endif

#ifdef __cplusplus
namespace amrex {
    /** Floating Point Type for Fields
     *
     * This is a floating point type used by default for fields such as
     * FArrayBox and MultiFab classes. This can be changed at compile time to
     * run double, single or experimental floating point precisions for field
     * computations and storage.
     */
    using Real = amrex_real;

    /** Floating Point Type for Particles
     *
     * This is a floating point type used by default for ParticleContainer and
     * related particle classes. This can be changed at compile time to run
     * double, single or experimental floating point precisions for particle
     * computations and storage.
     */
    using ParticleReal = amrex_particle_real;

/** Floating point literals for AMReX precision.
 *
 * This namespace can be used to write typed numeric literals in expressions.
 */
inline namespace literals {

    /** @{
      C++ user literals ``_rt`` &  ``_prt`` for short-hand notations

      Use this to properly add types to constant such as
      ```C++
      auto const mypi = 3.14_rt;
      auto const sphere_volume = 4_rt / 3_rt * pow(r, 3) * mypi;
      ```
    */
    constexpr Real
    operator"" _rt( long double x )
    {
        return Real( x );
    }

    constexpr Real
    operator"" _rt( unsigned long long int x )
    {
        return Real( x );
    }

    constexpr ParticleReal
    operator"" _prt( long double x )
    {
        return ParticleReal( x );
    }

    constexpr ParticleReal
    operator"" _prt( unsigned long long int x )
    {
        return ParticleReal( x );
    }
    /// @}
} // namespace literals
} // namespace amrex
#endif

#endif /* !BL_LANG_FORT */

#endif /*BL_REAL_H*/
