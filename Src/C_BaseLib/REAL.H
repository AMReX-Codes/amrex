
#ifndef BL_REAL_H
#define BL_REAL_H

/*
** $Id: REAL.H,v 1.5 2000-10-02 20:52:37 lijewski Exp $
*/

#if defined(BL_ARCH_CRAY)
#  ifdef BL_USE_DOUBLE
#    error BL_USE_DOUBLE not allowed on CRAY architectures
#  endif
#  define BL_USE_FLOAT 1
#endif

#ifdef BL_USE_FLOAT
#  undef BL_USE_DOUBLE
#  undef BL_USE_FLOAT
#  define BL_USE_FLOAT 1
#  if !defined(BL_LANG_FORT)

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

/*@ManDoc:
  REAL is a typedef specifying the precision of the floating-point
  calculations in C++ code.  It will be either `float' or `double'
  depending upon which of the macros BL\_USE\_FLOAT or
  BL\_USE\_DOUBLE, respectively, is defined during compilations.  For
  portability, you should write floating-point code in terms of this
  typedef, instead of using `float' or `double' directly.

  Note that exactly one of these macros must be defined
  when compiling any module that uses floating-point.  Also,
  on Cray architectures all floating-point computations are done in 64 bits.
*/
     typedef float REAL;

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
     typedef float Real;
#  else             /* LANG_FORT */
#    define REAL_T real*4
#  endif
#else               /* ! USE_FLOAT */
#  undef  BL_USE_FLOAT
#  undef  BL_USE_DOUBLE
#  define BL_USE_DOUBLE 1
#  if !defined(BL_LANG_FORT)
     typedef double REAL;
     typedef double Real;
#  else             /* LANG_FORT */

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
#    define REAL_T real*8
#  endif
#endif

#if !defined(BL_USE_FLOAT) && !defined(BL_USE_DOUBLE)
#error Exactly one of BL_USE_FLOAT or BL_USE_DOUBLE must be defined
#endif

#endif /*BL_REAL_H*/
