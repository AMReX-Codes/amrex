
#ifndef BL_ARRAYLIM_H
#define BL_ARRAYLIM_H

/*
** $Id: ArrayLim.H,v 1.7 2001-02-14 17:57:57 lijewski Exp $
*/

#if !defined(BL_LANG_FORT)

/*
** C stuff ...
*/

#ifdef BL_USE_ARLIM

#define ARLIM_P(x)  const int*
#define ARLIM(x)  x

#else

#if     BL_SPACEDIM==1
#define ARLIM_P(x)  const int&
#define ARLIM(x)  x[0]
#elif   BL_SPACEDIM==2
#define ARLIM_P(x)  const int&,const int&
#define ARLIM(x)  x[0],x[1]
#elif   BL_SPACEDIM==3
#define ARLIM_P(x)  const int&,const int&,const int&
#define ARLIM(x)  x[0],x[1],x[2]
#endif
#endif

#define BL_TO_FORTRAN(x) (x).dataPtr(), ARLIM((x).loVect()), ARLIM((x).hiVect())
#define BL_TO_FORTRAN_N(x,n) (x).dataPtr(n), ARLIM((x).loVect()), ARLIM((x).hiVect())

#else

/*
** Fortran stuff ...
*/

#include <REAL.H>

#ifdef BL_USE_ARLIM

#error "This section of code is deprecated: BL_USE_ARLIM is invalid"

#if (BL_SPACEDIM == 1)
#define DIMS(a) a##_lo, a##_hi
#define DIMDEC(a) a##_lo(1), a##_hi(1)
#define DIMV(a) a##_lo(1):a##_hi(1)
#define DIM1(a) a##_lo(1):a##_hi(1)
#define ARG_L1(a) a##_lo(1)
#define ARG_H1(a) a##_hi(1)
#endif

#if (BL_SPACEDIM == 2)
#define DIMS(a) a##_lo, a##_hi
#define DIMDEC(a) a##_lo(2), a##_hi(2)
#define DIMV(a) a##_lo(1):a##_hi(1), a##_lo(2):a##_hi(2)
#define DIM1(a) a##_lo(1):a##_hi(1)
#define DIM2(a) a##_lo(2):a##_hi(2)
#define ARG_L1(a) a##_lo(1)
#define ARG_L2(a) a##_lo(2)
#define ARG_H1(a) a##_hi(1)
#define ARG_H2(a) a##_hi(2)
#endif

#if (BL_SPACEDIM == 3)
#define DIMS(a) a##_lo, a##_hi
#define DIMDEC(a) a##_lo(3), a##_hi(3)
#define DIMV(a) a##_lo(1):a##_hi(1), a##_lo(2):a##_hi(2), a##_lo(3):a##_hi(3)
#define DIM1(a) a##_lo(1):a##_hi(1)
#define DIM2(a) a##_lo(2):a##_hi(2)
#define DIM3(a) a##_lo(3):a##_hi(3)
#define DIM12(a) a##_lo(1):a##_hi(1), a##_lo(2):a##_hi(2)
#define DIM23(a) a##_lo(2):a##_hi(2), a##_lo(3):a##_hi(3)
#define DIM13(a) a##_lo(1):a##_hi(1), a##_lo(3):a##_hi(3)
#define ARG_L1(a) a##_lo(1)
#define ARG_L2(a) a##_lo(2)
#define ARG_L3(a) a##_lo(3)
#define ARG_H1(a) a##_hi(1)
#define ARG_H2(a) a##_hi(2)
#define ARG_H3(a) a##_hi(3)
#endif

#else

#if __STDC__==1

#if (BL_SPACEDIM == 1)
#define DIMS(a) a##_l1, a##_h1
#define DIMDEC(a) a##_l1, a##_h1
#define DIMV(a) a##_l1:a##_h1
#define DIM1(a) a##_l1:a##_h1
#define ARG_L1(a) a##_l1
#define ARG_H1(a) a##_h1
#endif

#if (BL_SPACEDIM == 2)
#define DIMS(a) a##_l1, a##_l2 , a##_h1, a##_h2
#define DIMDEC(a) a##_l1, a##_l2 , a##_h1, a##_h2
#define DIMV(a) a##_l1:a##_h1, a##_l2:a##_h2
#define DIM1(a) a##_l1:a##_h1
#define DIM2(a) a##_l2:a##_h2
#define ARG_L1(a) a##_l1
#define ARG_L2(a) a##_l2
#define ARG_H1(a) a##_h1
#define ARG_H2(a) a##_h2
#endif

#if (BL_SPACEDIM == 3)
#define DIMS(a) a##_l1, a##_l2 , a##_l3, a##_h1, a##_h2, a##_h3
#define DIMDEC(a) a##_l1, a##_l2 , a##_l3, a##_h1, a##_h2, a##_h3
#define DIMV(a) a##_l1:a##_h1, a##_l2:a##_h2, a##_l3:a##_h3
#define DIM1(a) a##_l1:a##_h1
#define DIM2(a) a##_l2:a##_h2
#define DIM3(a) a##_l3:a##_h3
#define DIM12(a) a##_l1:a##_h1, a##_l2:a##_h2
#define DIM23(a) a##_l2:a##_h2, a##_l3:a##_h3
#define DIM13(a) a##_l1:a##_h1, a##_l3:a##_h3
#define ARG_L1(a) a##_l1
#define ARG_L2(a) a##_l2
#define ARG_L3(a) a##_l3
#define ARG_H1(a) a##_h1
#define ARG_H2(a) a##_h2
#define ARG_H3(a) a##_h3
#endif

#else

#if (BL_SPACEDIM == 1)
#define DIMS(a) a/**/_l1, a/**/_h1
#define DIMDEC(a) a/**/_l1, a/**/_h1
#define DIMV(a) a/**/_l1:a/**/_h1
#define DIM1(a) a/**/_l1:a/**/_h1
#define ARG_L1(a) a/**/_l1
#define ARG_H1(a) a/**/_h1
#endif

#if (BL_SPACEDIM == 2)
#define DIMS(a) a/**/_l1, a/**/_l2 , a/**/_h1, a/**/_h2
#define DIMDEC(a) a/**/_l1, a/**/_l2 , a/**/_h1, a/**/_h2
#define DIMV(a) a/**/_l1:a/**/_h1, a/**/_l2:a/**/_h2
#define DIM1(a) a/**/_l1:a/**/_h1
#define DIM2(a) a/**/_l2:a/**/_h2
#define ARG_L1(a) a/**/_l1
#define ARG_L2(a) a/**/_l2
#define ARG_H1(a) a/**/_h1
#define ARG_H2(a) a/**/_h2
#endif

#if (BL_SPACEDIM == 3)
#define DIMS(a) a/**/_l1, a/**/_l2 , a/**/_l3, a/**/_h1, a/**/_h2, a/**/_h3
#define DIMDEC(a) a/**/_l1, a/**/_l2 , a/**/_l3, a/**/_h1, a/**/_h2, a/**/_h3
#define DIMV(a) a/**/_l1:a/**/_h1, a/**/_l2:a/**/_h2, a/**/_l3:a/**/_h3
#define DIM1(a) a/**/_l1:a/**/_h1
#define DIM2(a) a/**/_l2:a/**/_h2
#define DIM3(a) a/**/_l3:a/**/_h3
#define DIM12(a) a/**/_l1:a/**/_h1, a/**/_l2:a/**/_h2
#define DIM23(a) a/**/_l2:a/**/_h2, a/**/_l3:a/**/_h3
#define DIM13(a) a/**/_l1:a/**/_h1, a/**/_l3:a/**/_h3
#define ARG_L1(a) a/**/_l1
#define ARG_L2(a) a/**/_l2
#define ARG_L3(a) a/**/_l3
#define ARG_H1(a) a/**/_h1
#define ARG_H2(a) a/**/_h2
#define ARG_H3(a) a/**/_h3
#endif

#endif /*__STDC__==1*/

#endif /*BL_USE_ARLIM*/

#endif /*!defined(BL_LANG_FORT)*/

#endif /*BL_ARRAYLIM_H*/
