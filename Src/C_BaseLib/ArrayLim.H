//BL_COPYRIGHT_NOTICE

#ifndef _ARRAYLIM_H_
#define _ARRAYLIM_H_

//
// $Id: ArrayLim.H,v 1.1 1997-09-12 18:00:00 lijewski Exp $
//

#if !defined(BL_LANG_FORT)

//
// The #ifdef ... #endif subterfuge here is to force doc++ to properly
// document BL_USE_ARLIM :-(
//
#ifdef BL_USE_ARLIM
#undef BL_USE_ARLIM
/*@ManDoc:
  The macro BL\_USE\_ARLIM indicates that the Fortran compiler allows array
  limits to be passed as integer vectors.  The default for Standard Fortran
  is to expand the integer vector into its respective components.  This
  macro is used to control the expansion of the macros ARLIM\_P(x) and
  ARLIM(x).
*/
#define BL_USE_ARLIM 1

/*@ManDoc:
  The macro ARLIM\_P(x) should be used in C++ code when prototyping a Fortran
  function in an extern "C" statement (or block) that takes a
  BL\_SPACEDIM-dimensional array and a BL\_SPACEDIM-dimensional integer vector
  indicating the upper and lower dimensions of the array.

  For example:

    extern "C" void f (Real* array, ARLIM\_P(ivec));

  When BL\_USE\_ARLIM is defined, this macro expands to

    const int *

  When BL\_USE\_ARLIM is not defined, this macro expands to one of

    const int \&

    const int \&, const int \&

    const int \&, const int \&, const int \&

  depending on the value of BL\_SPACEDIM.
*/
# define ARLIM\_P(x)  const int*

/*@ManDoc:
  The macro ARLIM(x) should be used in C++ code when passing an integer vector
  specifying the dimensions of an array, to a Fortran routine that has
  previously been prototyped with ARLIM\_P(x).

  For example:

     f(array, ARLIM(ivec));

  When BL\_USE\_ARLIM is defined, this macro expands to

    x

  When BL\_USE\_ARLIM is not defined, this macro expands to one of

    x[0]

    x[0],x[1]

    x[0],x[1],x[2]

  depending on the value of BL\_SPACEDIM.
*/
# define ARLIM(x)  x
#else
#   if   BL_SPACEDIM==1
#       define  ARLIM_P(x)  const int&
#       define  ARLIM(x)  x[0]
#   elif BL_SPACEDIM==2
#       define  ARLIM_P(x)  const int&,const int&
#       define  ARLIM(x)  x[0],x[1]
#   elif BL_SPACEDIM==3
#       define  ARLIM_P(x)  const int&,const int&,const int&
#       define  ARLIM(x)  x[0],x[1],x[2]
#   endif
#endif

#define BL_FAA_TYPE(u) REAL* u, ARLIM_P(u), ARLIM_P(u)

#define BL_BOXARG_TYPE(bx) const int*, const int*

#define BL_FAA(u) u.dataPtr(),ARLIM(u.loVect()),ARLIM(u.hiVect())

#define BL_FAAN(u,n) u.dataPtr(n),ARLIM(u.loVect()),ARLIM(u.hiVect())

#define BL_BOXARG(b) b.loVect(), b.hiVect()

#else /* Fortran stuff */

#include <REAL.H>

#if (BL_SPACEDIM == 1)
#define BL_FARG(u) u, u##loi1, u##hii1
#define BL_BBOUNDS(b) integer b##lo(1), b##hi(1)
#define BL_FBOUNDS(u) integer u##loi1, u##hii1
#define BL_FARRAY(u,n) REAL_T u(u##loi1:u##hii1, n)
#define BL_FARRAY1(u) REAL_T u(u##loi1:u##hii1)
#endif

#if (BL_SPACEDIM == 2)
#define BL_FARG(u) u, u##loi1, u##loi2, u##hii1, u##hii2
#define BL_BBOUNDS(b) integer b##lo(2), b##hi(2)
#define BL_FBOUNDS(u) integer u##loi1, u##loi2, u##hii1, u##hii2
#define BL_FARRAY(u,n) REAL_T u(u##loi1:u##hii1, u##loi2:u##hii2, n)
#define BL_FARRAY1(u) REAL_T u(u##loi1:u##hii1, u##loi2:u##hii2)
#endif

#if (BL_SPACEDIM == 3)
#define BL_FARG(u) u, u##loi1, u##loi2, u##loi3, u##hii1, u##hii2, u##hii3
#define BL_BBOUNDS(b) integer b##lo(3), b##hi(3)
#define BL_FBOUNDS(u) integer u##loi1, u##loi2, u##loi3, u##hii1, u##hii2, u##hii3
#define BL_FARRAY(u,n) REAL_T u(u##loi1:u##hii1, u##loi2:u##hii2, u##loi3:u##hii3, n)
#define BL_FARRAY1(u) REAL_T u(u##loi1:u##hii1, u##loi2:u##hii2, u##loi3:u##hii3)
#endif

#define BL_BARG(b) b##lo, b##hi

#endif

#endif /*_ARRAYLIM_H_*/
