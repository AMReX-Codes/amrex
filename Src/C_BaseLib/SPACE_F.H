/* BL_COPYRIGHT_NOTICE */

#ifndef BL_SPACE_F_H
#define BL_SPACE_F_H

/*
/* $Id: SPACE_F.H,v 1.3 1997-09-18 18:03:34 car Exp $
/*

#include <SPACE.H>

#ifdef BL_LANG_FORT
#  if (BL_SPACEDIM == 1)

/*@ManDoc:
  The DIM\_V() macro can be used to enhance the portability of Fortran
  code.  It can be used whenever the BL\_LANG\_FORT preprocessor macro
  is defined.  It is used to specify the dimensionality of Fortran
  arrays, corresponding to BL\_SPACEDIM.  `lo' and `hi' are
  BL\_SPACEDIM-dimensional integer arrays that specify the low and
  high range of the array being specified.

  More explicitly, the macro evaluates to one of


     lo(1):hi(1)

     lo(1):hi(1),lo(2):hi(2)

     lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)

  corresponding to BL\_SPACEDIM values of 1, 2, or 3.
*/
#    define DIM_V(lo,hi) lo(1):hi(1)
#  elif (BL_SPACEDIM == 2)
#    define DIM_V(lo,hi) lo(1):hi(1),lo(2):hi(2)
#  elif (BL_SPACEDIM == 3)
#    define DIM_V(lo,hi) lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)
#  endif
#endif

#endif /*BL_SPACE_F_H*/
