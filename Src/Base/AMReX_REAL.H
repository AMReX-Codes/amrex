
#ifndef AMREX_REAL_H
#define AMREX_REAL_H

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
#else
typedef double amrex_real;
#endif

#ifdef __cplusplus
namespace amrex {
    using Real = amrex_real;
}
#endif

#else

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

#ifdef NAGFOR
#ifndef CONCAT
#define IDENTITY(a) a
#define CONCAT(a,b) IDENTITY(a)IDENTITY(b)
#endif
#endif

#ifdef BL_USE_FLOAT
#    define REAL_T REAL
#if __STDC__==1 || defined(__INTEL_COMPILER)
#        define BL_REAL(a) a##E0
#        define BL_REAL_E(a,b) a##E##b
#elif NAGFOR
#        define BL_REAL(a) CONCAT(a,E0)
#        define BL_REAL_E(a,b) CONCAT(,CONCAT(E,b))
#else
#        define BL_REAL(a) a/**/E0
#        define BL_REAL_E(a,b) a/**/E/**/b
#endif
#else
#    define REAL_T DOUBLE PRECISION
#if __STDC__==1 || defined(__INTEL_COMPILER)
#        define BL_REAL(a) a##D0
#        define BL_REAL_E(a,b) a##D##b
#elif NAGFOR
#        define BL_REAL(a) CONCAT(a,D0)
#        define BL_REAL_E(a,b) CONCAT(,CONCAT(D,b))
#else
#        define BL_REAL(a) a/**/D0
#        define BL_REAL_E(a,b) a/**/D/**/b
#endif
#endif

#endif /* !BL_LANG_FORT */

#endif /*BL_REAL_H*/
