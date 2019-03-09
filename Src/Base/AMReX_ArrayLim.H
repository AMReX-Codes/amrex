#ifndef AMREX_ARRAYLIM_H_
#define AMREX_ARRAYLIM_H_

#include <AMReX_REAL.H>

#ifdef __cplusplus

#include <AMReX_Array.H>

/*
** C++ stuff ...
*/

#ifndef AMREX_XSDK

#if     AMREX_SPACEDIM==1
#define ARLIM_P(x)  const int&
#define ARLIM(x)  (x)[0]
#define ARLIM_3D(x) amrex::GpuArray<int,3>{(x)[0], 0, 0}.data()
#define ARLIM_ANYD(x) amrex::GpuArray<int,3>{(x)[0], 0, 0}.data()
#define ZFILL(x) amrex::GpuArray<amrex::Real,3>{(x)[0], 0., 0.}.data()
#ifdef AMREX_USE_GPU_PRAGMA
#define AMREX_REAL_ANYD(x) (x)[0], 0., 0.
#else
#define AMREX_REAL_ANYD(x) AMREX_ZFILL(x)
#endif
#elif   AMREX_SPACEDIM==2
#define ARLIM_P(x)  const int&,const int&
#define ARLIM(x)  (x)[0],(x)[1]
#define ARLIM_3D(x) amrex::GpuArray<int,3>{(x)[0], (x)[1], 0}.data()
#define ARLIM_ANYD(x) amrex::GpuArray<int,3>{(x)[0], (x)[1], 0}.data()
#define ZFILL(x) amrex::GpuArray<amrex::Real,3>{(x)[0], (x)[1], 0.}.data()
#ifdef AMREX_USE_GPU_PRAGMA
#define AMREX_REAL_ANYD(x) (x)[0], (x)[1], 0.
#else
#define AMREX_REAL_ANYD(x) AMREX_ZFILL(x)
#endif
#elif   AMREX_SPACEDIM==3
#define ARLIM_P(x)  const int&,const int&,const int&
#define ARLIM(x)  (x)[0],(x)[1],(x)[2]
#define ARLIM_3D(x) x
#define ARLIM_ANYD(x) x
#define ZFILL(x) x
#ifdef AMREX_USE_GPU_PRAGMA
#define AMREX_REAL_ANYD(x) (x)[0], (x)[1], (x)[2]
#else
#define AMREX_REAL_ANYD(x) AMREX_ZFILL(x)
#endif
#endif

#endif /* ndef AMREX_XSDK */

#if     AMREX_SPACEDIM==1
#define AMREX_ARLIM_P(x)  const int&
#define AMREX_ARLIM(x)  (x)[0]
#define AMREX_BCREC_3D(x) amrex::GpuArray<int,6>{(x)[0], 0, 0, (x)[1], 0, 0}.data()
#define AMREX_BCREC_ANYD(x) amrex::GpuArray<int,6>{(x)[0], 0, 0, (x)[1], 0, 0}.data()
#define AMREX_ARLIM_3D(x) amrex::GpuArray<int,3>{(x)[0], 0, 0}.data()
#define AMREX_ARLIM_ANYD(x) amrex::GpuArray<int,3>{(x)[0], 0, 0}.data()
#define AMREX_ZFILL(x) amrex::GpuArray<amrex::Real,3>{(x)[0], 0., 0.}.data()
#ifdef AMREX_USE_GPU_PRAGMA
#define AMREX_ARLIM_ARG(x) x[0], 0, 0
#else
#define AMREX_ARLIM_ARG(x) AMREX_ARLIM_3D(x)
#endif
#elif   AMREX_SPACEDIM==2
#define AMREX_ARLIM_P(x)  const int&,const int&
#define AMREX_ARLIM(x)  (x)[0],(x)[1]
#define AMREX_BCREC_3D(x) amrex::GpuArray<int,6>{(x)[0], (x)[1], 0, (x)[2], (x)[3], 0}.data()
#define AMREX_BCREC_ANYD(x) amrex::GpuArray<int,6>{(x)[0], (x)[1], 0, (x)[2], (x)[3], 0}.data()
#define AMREX_ARLIM_3D(x) amrex::GpuArray<int,3>{(x)[0], (x)[1], 0}.data()
#define AMREX_ARLIM_ANYD(x) amrex::GpuArray<int,3>{(x)[0], (x)[1], 0}.data()
#define AMREX_ZFILL(x) amrex::GpuArray<amrex::Real,3>{(x)[0], (x)[1], 0.}.data()
#ifdef AMREX_USE_GPU_PRAGMA
#define AMREX_ARLIM_ARG(x) x[0], x[1], 0
#else
#define AMREX_ARLIM_ARG(x) AMREX_ARLIM_3D(x)
#endif
#elif   AMREX_SPACEDIM==3
#define AMREX_ARLIM_P(x)  const int&,const int&,const int&
#define AMREX_ARLIM(x)  (x)[0],(x)[1],(x)[2]
#define AMREX_ARLIM_3D(x) x
#define AMREX_ARLIM_ANYD(x) x
#define AMREX_ZFILL(x) x
#define AMREX_BCREC_3D(x) x
#define AMREX_BCREC_ANYD(x) x
#ifdef AMREX_USE_GPU_PRAGMA
#define AMREX_ARLIM_ARG(x) x[0], x[1], x[2]
#else
#define AMREX_ARLIM_ARG(x) AMREX_ARLIM_3D(x)
#endif
#endif

#ifndef AMREX_USE_GPU_PRAGMA
#define AMREX_INT_ANYD(x) AMREX_ARLIM_ANYD(x)
#else
#define AMREX_INT_ANYD(x) AMREX_ARLIM_ARG(x)
#endif

#define BL_TO_FORTRAN(x) (x).dataPtr(), AMREX_ARLIM((x).loVect()), AMREX_ARLIM((x).hiVect())
#define BL_TO_FORTRAN_FAB(x) (x).dataPtr(), AMREX_ARLIM_3D((x).loVect()), AMREX_ARLIM_3D((x).hiVect()), (x).nCompPtr()
#define BL_TO_FORTRAN_BOX(x) AMREX_ARLIM_3D((x).loVect()), AMREX_ARLIM_3D((x).hiVect())
#define BL_TO_FORTRAN_N(x,n) (x).dataPtr(n), AMREX_ARLIM((x).loVect()), AMREX_ARLIM((x).hiVect())

#define BL_TO_FORTRAN_3D(x) (x).dataPtr(), AMREX_ARLIM_3D((x).loVect()), AMREX_ARLIM_3D((x).hiVect())
#define BL_TO_FORTRAN_ANYD(x) (x).dataPtr(), AMREX_ARLIM_ANYD((x).loVect()), AMREX_ARLIM_ANYD((x).hiVect())
#define BL_TO_FORTRAN_N_3D(x,n) (x).dataPtr(n), AMREX_ARLIM_3D((x).loVect()), AMREX_ARLIM_3D((x).hiVect())
#define BL_TO_FORTRAN_N_ANYD(x,n) (x).dataPtr(n), AMREX_ARLIM_ANYD((x).loVect()), AMREX_ARLIM_ANYD((x).hiVect())

// Note that these GPU-friendly versions only cover the dimension-agnostic
// versions at this time. These macros are not intended for use by the developer
// but rather are inserted by the GPU code generation script.

#ifdef AMREX_USE_GPU_PRAGMA
#define BL_TO_FORTRAN_GPU_FAB(x) (x).dataPtr(), AMREX_INT_ANYD((x).loVect()), AMREX_INT_ANYD((x).hiVect()), (x).nCompPtr()
#define BL_TO_FORTRAN_GPU_BOX(x) AMREX_INT_ANYD((x).loVect()), AMREX_INT_ANYD((x).hiVect())
#define BL_TO_FORTRAN_GPU_3D(x) (x).dataPtr(), AMREX_INT_ANYD((x).loVect()), AMREX_INT_ANYD((x).hiVect())
#define BL_TO_FORTRAN_GPU_ANYD(x) (x).dataPtr(), AMREX_INT_ANYD((x).loVect()), AMREX_INT_ANYD((x).hiVect())
#define BL_TO_FORTRAN_GPU_N_3D(x,n) (x).dataPtr(n), AMREX_INT_ANYD((x).loVect()), AMREX_INT_ANYD((x).hiVect())
#define BL_TO_FORTRAN_GPU_N_ANYD(x,n) (x).dataPtr(n), AMREX_INT_ANYD((x).loVect()), AMREX_INT_ANYD((x).hiVect())
#endif

#else

#if !defined(BL_LANG_FORT)

#ifndef AMREX_XSDK

/* C stuff */
#if AMREX_SPACEDIM == 1
#define ARLIM_P(x) const int*
#elif AMREX_SPACEDIM == 2
#define ARLIM_P(x) const int*, const int*
#else
#define ARLIM_P(x) const int*, const int*, const int*
#endif

#endif /* ndef AMREX_XSDK */

#if AMREX_SPACEDIM == 1
#define AMREX_ARLIM_P(x) const int*
#elif AMREX_SPACEDIM == 2
#define AMREX_ARLIM_P(x) const int*, const int*
#else
#define AMREX_ARLIM_P(x) const int*, const int*, const int*
#endif

#else

#ifndef AMREX_XSDK

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

#endif /* ndef AMREX_XSDK */

#endif /*else of !defined(BL_LANG_FORT)*/

#endif /*else of __cplusplus*/

#endif /*BL_ARRAYLIM_H*/
