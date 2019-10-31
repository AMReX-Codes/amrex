#ifndef AMREX_REAL_H_
#define AMREX_REAL_H_

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
typedef float amrex_real;
// We need to define these to get around a CUDA 9.2 bug that breaks std::numeric_limits
#define AMREX_REAL_MIN     FLT_MIN
#define AMREX_REAL_MAX     FLT_MAX
#define AMREX_REAL_LOWEST -FLT_MAX
#else
typedef double amrex_real;
#define AMREX_REAL_MIN     DBL_MIN
#define AMREX_REAL_MAX     DBL_MAX
#define AMREX_REAL_LOWEST -DBL_MAX
#endif

#ifdef BL_SINGLE_PRECISION_PARTICLES
typedef float amrex_particle_real;
#else
typedef double amrex_particle_real;
#endif

#ifdef __cplusplus
namespace amrex {
    using Real = amrex_real;
    using ParticleReal = amrex_particle_real;

    /** @{
      C++ user literals ``_rt`` &  ``_prt`` for short-hand notations

      Use this to properly add types to constant such as
      ```C++
      auto const mypi = 3.14_rt;
      auto const sphere_volume = 4_rt / 3_rt * pow(r, 3) * mypi;
      ```
    */
    constexpr Real
    operator""_rt( long double x )
    {
        return Real( x );
    }

    constexpr Real
    operator""_rt( unsigned long long int x )
    {
        return Real( x );
    }

    constexpr ParticleReal
    operator""_prt( long double x )
    {
        return ParticleReal( x );
    }

    constexpr ParticleReal
    operator""_prt( unsigned long long int x )
    {
        return ParticleReal( x );
    }
    /// @}
} // namespace amrex
#endif

#else

#ifndef AMREX_XSDK

/*
  The REAL_T macro specifies the precision of the floating-point
  calculations in Fortran code.  It will be either real*4 or
  real*8 depending upon which of the symbols BL_USE_FLOAT or
  BL_USE_DOUBLE, respectively, is defined during compilations.  For
  portability, you should write floating-point code in terms of this
  macro, instead of using real*4 or real*8 directly.

  Note that exactly one of these macros must be defined
  when compiling any module that uses floating-point.
*/

#ifdef BL_USE_FLOAT
#    define REAL_T REAL
#if __STDC__==1 || defined(__INTEL_COMPILER)
#        define BL_REAL(a) a##E0
#        define BL_REAL_E(a,b) a##E##b
#else
#        define BL_REAL(a) a/**/E0
#        define BL_REAL_E(a,b) a/**/E/**/b
#endif
#else
#    define REAL_T DOUBLE PRECISION
#if __STDC__==1 || defined(__INTEL_COMPILER)
#        define BL_REAL(a) a##D0
#        define BL_REAL_E(a,b) a##D##b
#else
#        define BL_REAL(a) a/**/D0
#        define BL_REAL_E(a,b) a/**/D/**/b
#endif
#endif

#endif /* ndef AMREX_XSDK */

#endif /* !BL_LANG_FORT */

#endif /*BL_REAL_H*/
