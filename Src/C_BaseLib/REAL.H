#ifndef BL_REAL_H
#define BL_REAL_H

/*
 * $Id: REAL.H,v 1.11 2001-08-16 17:57:43 lijewski Exp $
 */

#ifdef BL_USE_FLOAT
#  undef BL_USE_DOUBLE
#  undef BL_USE_FLOAT
/*@ManDoc:
  The macro BL\_USE\_FLOAT indicates that C++ floating-point calculations should
  use "float" variables and Fortran floating-point calculations should use
  "real*4" variables.  One of BL\_USE\_FLOAT or BL\_USE\_DOUBLE must always be
  defined when compiling and using BoxLib.
*/
#  define BL_USE_FLOAT 1
#else
#  undef  BL_USE_FLOAT
#  undef  BL_USE_DOUBLE
/*@ManDoc:
  The macro BL\_USE\_DOUBLE indicates that C++ floating-point calculations
  should use "double" variables and Fortran floating-point calculations
  should use "real*8" variables.  One of BL\_USE\_FLOAT or BL\_USE\_DOUBLE must
  always be defined when compiling and using BoxLib.  This macro is not
  allowed on Cray architectures; i.e. only BL\_USE\_FLOAT is allowed on Cray
  architectures.
*/
#  define BL_USE_DOUBLE 1
#endif

#if !defined(BL_LANG_FORT)

/*@ManDoc: 
  Real is a typedef specifying the precision of the floating-point
  calculations in C++ code.  It will be either `float' or `double'
  depending upon which of the macros BL\_USE\_FLOAT or
  BL\_USE\_DOUBLE, respectively, is defined during compilations.  For
  portability, you should write floating-point code in terms of this
  typedef, instead of using `float' or `double' directly.

  Note that exactly one of these macros must be defined
  when compiling any module that uses floating-point.  Also,
  on Cray architectures all floating-point computations are done in 64 bits.
*/
#ifdef BL_USE_FLOAT
    typedef float Real;
#else
    typedef double Real;
#endif

#else

/*@ManDoc:
  The REAL\_T macro specifies the precision of the floating-point
  calculations in Fortran code.  It will be either `real*4' or
  `real*8' depending upon which of the symbols BL\_USE\_FLOAT or
  BL\_USE\_DOUBLE, respectively, is defined during compilations.  For
  portability, you should write floating-point code in terms of this
  macro, instead of using `real*4' or `real*8' directly.

  Note that exactly one of these macros must be defined
  when compiling any module that uses floating-point.  Also,
  on Cray architectures all floating-point computations are done in 64 bits.
*/

#ifdef BL_USE_FLOAT
#    define REAL_T real*4
#if __STDC__==1
#        define BL_REAL(a) a##E0
#        define BL_REAL_E(a,b) a##E##b
#else
#        define BL_REAL(a) a/**/E0
#        define BL_REAL_E(a,b) a/**/E/**/b
#endif
#else
#    define REAL_T real*8
#if __STDC__==1
#        define BL_REAL(a) a##D0
#        define BL_REAL_E(a,b) a##D##b
#else
#        define BL_REAL(a) a/**/D0
#        define BL_REAL_E(a,b) a/**/D/**/b
#endif
#endif

#endif /* !BL_LANG_FORT */

#if defined(BL_ARCH_CRAY)
#  if !defined(BL_USE_FLOAT)
#    error BL_USE_FLOAT required on CRAY architectures
#  endif
#endif

#endif /*BL_REAL_H*/
