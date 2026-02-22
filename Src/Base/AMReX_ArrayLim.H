#ifndef AMREX_ARRAYLIM_H_
#define AMREX_ARRAYLIM_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_ArrayLim.H
 * \brief Macros for passing AMReX box/array bounds to legacy Fortran kernels.
 */

#include <AMReX_REAL.H>

#ifdef __cplusplus

#include <AMReX_Array.H>

#if     AMREX_SPACEDIM==1
#define AMREX_ARLIM_P(x)  const int&
#define AMREX_ARLIM(x)  (x)[0]
#define AMREX_BCREC_3D(x)   amrex::GpuArray<int,6>{(x)[0], 0, 0, (x)[1], 0, 0}.data()
#define AMREX_BCREC_ANYD(x) amrex::GpuArray<int,6>{(x)[0], 0, 0, (x)[1], 0, 0}.data()
#define AMREX_ARLIM_3D(x)   amrex::GpuArray<int,3>{(x)[0], 0, 0}.data()
#define AMREX_ARLIM_ANYD(x) amrex::GpuArray<int,3>{(x)[0], 0, 0}.data()
#define AMREX_ZFILL(x) amrex::GpuArray<amrex::Real,3>{(x)[0], 0., 0.}.data()
#define AMREX_ARLIM_ARG(x) AMREX_ARLIM_3D(x)
#elif   AMREX_SPACEDIM==2
#define AMREX_ARLIM_P(x)  const int&,const int&
#define AMREX_ARLIM(x)  (x)[0],(x)[1]
#define AMREX_BCREC_3D(x)   amrex::GpuArray<int,6>{(x)[0], (x)[1], 0, (x)[2], (x)[3], 0}.data()
#define AMREX_BCREC_ANYD(x) amrex::GpuArray<int,6>{(x)[0], (x)[1], 0, (x)[2], (x)[3], 0}.data()
#define AMREX_ARLIM_3D(x)   amrex::GpuArray<int,3>{(x)[0], (x)[1], 0}.data()
#define AMREX_ARLIM_ANYD(x) amrex::GpuArray<int,3>{(x)[0], (x)[1], 0}.data()
#define AMREX_ZFILL(x) amrex::GpuArray<amrex::Real,3>{(x)[0], (x)[1], 0.}.data()
#define AMREX_ARLIM_ARG(x) AMREX_ARLIM_3D(x)
#elif   AMREX_SPACEDIM==3
#define AMREX_ARLIM_P(x)  const int&,const int&,const int&
#define AMREX_ARLIM(x)  (x)[0],(x)[1],(x)[2]
#define AMREX_ARLIM_3D(x) x
#define AMREX_ARLIM_ANYD(x) x
#define AMREX_ZFILL(x) x
#define AMREX_BCREC_3D(x) x
#define AMREX_BCREC_ANYD(x) x
#define AMREX_ARLIM_ARG(x) AMREX_ARLIM_3D(x)
#endif

#define AMREX_INT_ANYD(x) AMREX_ARLIM_ANYD(x)

/**
 * \def BL_TO_FORTRAN(x)
 * \brief Expand a BaseFab-like object into `(data, lo, hi)` arguments in the active dimensionality.
 * \param x Fab-like object providing `dataPtr()`, `loVect()`, and `hiVect()`.
 */
#define BL_TO_FORTRAN(x)      (x).dataPtr(), AMREX_ARLIM     ((x).loVect()), AMREX_ARLIM     ((x).hiVect())
/**
 * \def BL_TO_FORTRAN_3D(x)
 * \brief Expand a BaseFab-like object into explicit 3D `(data, lo, hi)` arguments.
 * \param x Fab-like object providing `dataPtr()`, `loVect()`, and `hiVect()`.
 */
#define BL_TO_FORTRAN_3D(x)   (x).dataPtr(), AMREX_ARLIM_3D  ((x).loVect()), AMREX_ARLIM_3D  ((x).hiVect())
/**
 * \def BL_TO_FORTRAN_ANYD(x)
 * \brief Expand a BaseFab-like object into ANYD `(data, lo, hi)` arguments padded to 3D.
 * \param x Fab-like object providing `dataPtr()`, `loVect()`, and `hiVect()`.
 */
#define BL_TO_FORTRAN_ANYD(x) (x).dataPtr(), AMREX_ARLIM_ANYD((x).loVect()), AMREX_ARLIM_ANYD((x).hiVect())
/**
 * \def BL_TO_FORTRAN_FAB(x)
 * \brief Expand a BaseFab-like object into `(data, lo, hi, ncomp)` arguments for Fortran calls.
 * \param x Fab-like object providing `dataPtr()`, bounds, and `nCompPtr()`.
 */
#define BL_TO_FORTRAN_FAB(x)  (x).dataPtr(), AMREX_ARLIM_3D  ((x).loVect()), AMREX_ARLIM_3D  ((x).hiVect()), (x).nCompPtr()

/**
 * \def BL_TO_FORTRAN_N(x,n)
 * \brief As BL_TO_FORTRAN but passing the data pointer for component \p n.
 * \param x Fab-like object.
 * \param n Component index passed to `dataPtr(n)`.
 */
#define BL_TO_FORTRAN_N(x,n)      (x).dataPtr(n), AMREX_ARLIM     ((x).loVect()), AMREX_ARLIM     ((x).hiVect())
/// \copydoc BL_TO_FORTRAN_3D
#define BL_TO_FORTRAN_N_3D(x,n)   (x).dataPtr(n), AMREX_ARLIM_3D  ((x).loVect()), AMREX_ARLIM_3D  ((x).hiVect())
/// \copydoc BL_TO_FORTRAN_ANYD
#define BL_TO_FORTRAN_N_ANYD(x,n) (x).dataPtr(n), AMREX_ARLIM_ANYD((x).loVect()), AMREX_ARLIM_ANYD((x).hiVect())

/**
 * \def BL_TO_FORTRAN_BOX(x)
 * \brief Expand a Box into `(lo, hi)` arrays padded to 3D for Fortran kernels.
 * \param x Box supplying `loVect()` and `hiVect()`.
 */
#define BL_TO_FORTRAN_BOX(x) AMREX_ARLIM_3D((x).loVect()), AMREX_ARLIM_3D((x).hiVect())

#endif /* __cplusplus */

#endif /* AMREX_ARRAYLIM_H_ */
