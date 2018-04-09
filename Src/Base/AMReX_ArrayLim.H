
#ifndef BL_ARRAYLIM_H
#define BL_ARRAYLIM_H

#include <AMReX_REAL.H>

#ifdef __cplusplus

#include <array>

/*
** C++ stuff ...
*/

#if     AMREX_SPACEDIM==1
#define ARLIM_P(x)  const int&
#define ARLIM(x)  (x)[0]
#define ARLIM_3D(x) std::array<int,3>{(x)[0], 0, 0}.data()
#define ARLIM_ANYD(x) std::array<int,3>{(x)[0], 0, 0}.data()
#define ZFILL(x) std::array<amrex::Real,3>{(x)[0], 0., 0.}.data()
#define BCREC_3D(x) std::array<int,6>{(x)[0], 0, 0, (x)[1], 0, 0}.data()
#define BCREC_ANYD(x) std::array<int,6>{(x)[0], 0, 0, (x)[1], 0, 0}.data()
#elif   AMREX_SPACEDIM==2
#define ARLIM_P(x)  const int&,const int&
#define ARLIM(x)  (x)[0],(x)[1]
#define ARLIM_3D(x) std::array<int,3>{(x)[0], (x)[1], 0}.data()
#define ARLIM_ANYD(x) std::array<int,3>{(x)[0], (x)[1], 0}.data()
#define ZFILL(x) std::array<amrex::Real,3>{(x)[0], (x)[1], 0.}.data()
#define BCREC_3D(x) std::array<int,6>{(x)[0], (x)[1], 0, (x)[2], (x)[3], 0}.data()
#define BCREC_ANYD(x) std::array<int,6>{(x)[0], (x)[1], 0, (x)[2], (x)[3], 0}.data()
#elif   AMREX_SPACEDIM==3
#define ARLIM_P(x)  const int&,const int&,const int&
#define ARLIM(x)  (x)[0],(x)[1],(x)[2]
#define ARLIM_3D(x) x
#define ARLIM_ANYD(x) x
#define ZFILL(x) x
#define BCREC_3D(x) x
#define BCREC_ANYD(x) x
#endif

#define BL_TO_FORTRAN(x) (x).dataPtr(), ARLIM((x).loVect()), ARLIM((x).hiVect())
#define BL_TO_FORTRAN_FAB(x) (x).dataPtr(), ARLIM_3D((x).loVect()), ARLIM_3D((x).hiVect()), (x).nCompPtr()
#define BL_TO_FORTRAN_BOX(x) ARLIM_3D((x).loVect()), ARLIM_3D((x).hiVect())
#define BL_TO_FORTRAN_N(x,n) (x).dataPtr(n), ARLIM((x).loVect()), ARLIM((x).hiVect())

#define BL_TO_FORTRAN_3D(x) (x).dataPtr(), ARLIM_3D((x).loVect()), ARLIM_3D((x).hiVect())
#define BL_TO_FORTRAN_ANYD(x) (x).dataPtr(), ARLIM_ANYD((x).loVect()), ARLIM_ANYD((x).hiVect())
#define BL_TO_FORTRAN_N_3D(x,n) (x).dataPtr(n), ARLIM_3D((x).loVect()), ARLIM_3D((x).hiVect())
#define BL_TO_FORTRAN_N_ANYD(x,n) (x).dataPtr(n), ARLIM_ANYD((x).loVect()), ARLIM_ANYD((x).hiVect())

#else

#if !defined(BL_LANG_FORT)

/* C stuff */
#if AMREX_SPACEDIM == 1
#define ARLIM_P(x) const int*
#elif AMREX_SPACEDIM == 2
#define ARLIM_P(x) const int*, const int*
#else
#define ARLIM_P(x) const int*, const int*, const int*
#endif

#else

/*
** Fortran stuff ...
*/

#if __STDC__==1 || defined(__INTEL_COMPILER)

#if (AMREX_SPACEDIM == 1)
#define DIMS(a) a##_l1, a##_h1
#define DIMDEC(a) a##_l1, a##_h1
#define DIMV(a) a##_l1:a##_h1
#define DIM1(a) a##_l1:a##_h1
#define ARG_L1(a) a##_l1
#define ARG_H1(a) a##_h1
#define DIMARG(a) a##(1)
#endif

#if (AMREX_SPACEDIM == 2)
#define DIMS(a) a##_l1, a##_l2, a##_h1, a##_h2
#define DIMDEC(a) a##_l1, a##_l2, a##_h1, a##_h2
#define DIMV(a) a##_l1:a##_h1, a##_l2:a##_h2
#define DIM1(a) a##_l1:a##_h1
#define DIM2(a) a##_l2:a##_h2
#define ARG_L1(a) a##_l1
#define ARG_L2(a) a##_l2
#define ARG_H1(a) a##_h1
#define ARG_H2(a) a##_h2
#define DIMARG(a) a##(1),a##(2)
#endif

#if (AMREX_SPACEDIM == 3)
#define DIMS(a) a##_l1, a##_l2, a##_l3, a##_h1, a##_h2, a##_h3
#define DIMDEC(a) a##_l1, a##_l2, a##_l3, a##_h1, a##_h2, a##_h3
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
#define DIMARG(a) a##(1),a##(2),a##(3)
#endif

#elif NAGFOR

#ifndef CONCAT
#define IDENTITY(a) a
#define CONCAT(a,b) IDENTITY(a)IDENTITY(b)
#endif

#if (AMREX_SPACEDIM == 1)
#define DIMS(a) CONCAT(a,_l1), CONCAT(a,_h1)
#define DIMDEC(a) CONCAT(a,_l1), CONCAT(a,_h1)
#define DIMV(a) CONCAT(a,_l1):CONCAT(a,_h1)
#define DIM1(a) CONCAT(a,_l1):CONCAT(a,_h1)
#define ARG_L1(a) CONCAT(a,_l1)
#define ARG_H1(a) CONCAT(a,_h1)
#define DIMARG(a) CONCAT(a,(1))
#endif

#if (AMREX_SPACEDIM == 2)
#define DIMS(a) CONCAT(a,_l1), CONCAT(a,_l2), CONCAT(a,_h1), CONCAT(a,_h2)
#define DIMDEC(a) CONCAT(a,_l1), CONCAT(a,_l2), CONCAT(a,_h1), CONCAT(a,_h2)
#define DIMV(a) CONCAT(a,_l1):CONCAT(a,_h1), CONCAT(a,_l2):CONCAT(a,_h2)
#define DIM1(a) CONCAT(a,_l1):CONCAT(a,_h1)
#define DIM2(a) CONCAT(a,_l2):CONCAT(a,_h2)
#define ARG_L1(a) CONCAT(a,_l1)
#define ARG_L2(a) CONCAT(a,_l2)
#define ARG_H1(a) CONCAT(a,_h1)
#define ARG_H2(a) CONCAT(a,_h2)
#define DIMARG(a) CONCAT(a,(1)),CONCAT(a,(2))
#endif

#if (AMREX_SPACEDIM == 3)
#define DIMS(a) CONCAT(a,_l1), CONCAT(a,_l2), CONCAT(a,_l3), CONCAT(a,_h1), CONCAT(a,_h2), CONCAT(a,_h3)
#define DIMDEC(a) CONCAT(a,_l1), CONCAT(a,_l2), CONCAT(a,_l3), CONCAT(a,_h1), CONCAT(a,_h2), CONCAT(a,_h3)
#define DIMV(a) CONCAT(a,_l1):CONCAT(a,_h1), CONCAT(a,_l2):CONCAT(a,_h2), CONCAT(a,_l3):CONCAT(a,_h3)
#define DIM1(a) CONCAT(a,_l1):CONCAT(a,_h1)
#define DIM2(a) CONCAT(a,_l2):CONCAT(a,_h2)
#define DIM3(a) CONCAT(a,_l3):CONCAT(a,_h3)
#define DIM12(a) CONCAT(a,_l1):CONCAT(a,_h1), CONCAT(a,_l2):CONCAT(a,_h2)
#define DIM23(a) CONCAT(a,_l2):CONCAT(a,_h2), CONCAT(a,_l3):CONCAT(a,_h3)
#define DIM13(a) CONCAT(a,_l1):CONCAT(a,_h1), CONCAT(a,_l3):CONCAT(a,_h3)
#define ARG_L1(a) CONCAT(a,_l1)
#define ARG_L2(a) CONCAT(a,_l2)
#define ARG_L3(a) CONCAT(a,_l3)
#define ARG_H1(a) CONCAT(a,_h1)
#define ARG_H2(a) CONCAT(a,_h2)
#define ARG_H3(a) CONCAT(a,_h3)
#define DIMARG(a) CONCAT(a,(1)),CONCAT(a,(2)),CONCAT(a,(3))
#endif

#else

#if (AMREX_SPACEDIM == 1)
#define DIMS(a) a/**/_l1, a/**/_h1
#define DIMDEC(a) a/**/_l1, a/**/_h1
#define DIMV(a) a/**/_l1:a/**/_h1
#define DIM1(a) a/**/_l1:a/**/_h1
#define ARG_L1(a) a/**/_l1
#define ARG_H1(a) a/**/_h1
#define DIMARG(a) a/**/(1)
#endif

#if (AMREX_SPACEDIM == 2)
#define DIMS(a) a/**/_l1, a/**/_l2, a/**/_h1, a/**/_h2
#define DIMDEC(a) a/**/_l1, a/**/_l2, a/**/_h1, a/**/_h2
#define DIMV(a) a/**/_l1:a/**/_h1, a/**/_l2:a/**/_h2
#define DIM1(a) a/**/_l1:a/**/_h1
#define DIM2(a) a/**/_l2:a/**/_h2
#define ARG_L1(a) a/**/_l1
#define ARG_L2(a) a/**/_l2
#define ARG_H1(a) a/**/_h1
#define ARG_H2(a) a/**/_h2
#define DIMARG(a) a/**/(1),a/**/(2)
#endif

#if (AMREX_SPACEDIM == 3)
#define DIMS(a) a/**/_l1, a/**/_l2, a/**/_l3, a/**/_h1, a/**/_h2, a/**/_h3
#define DIMDEC(a) a/**/_l1, a/**/_l2, a/**/_l3, a/**/_h1, a/**/_h2, a/**/_h3
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
#define DIMARG(a) a/**/(1),a/**/(2),a/**/(3)
#endif

#endif /*__STDC__==1*/

#endif /*else of !defined(BL_LANG_FORT)*/

#endif /*else of __cplusplus*/

#endif /*BL_ARRAYLIM_H*/
