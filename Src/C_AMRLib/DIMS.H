#ifndef _DIMS_H_
#define _DIMS_H_ 

#ifdef BL_USE_ARLIM

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

#endif

#endif
