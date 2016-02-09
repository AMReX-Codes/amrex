#ifndef _ABec2_F_H_
#define _ABec2_F_H_

#include <BLFort.H>

BL_FORT_PROC_DECL(AB2_BNDRYRLX, ab2_bndryrlx)
    ( const int lo[], const int hi[],
      BL_FORT_FAB_ARG(f),
      const BL_FORT_IFAB_ARG(m),
      const int* cdir, const int* bct,
      const Real* bcl, const int* maxorder, const Real dx[] );

BL_FORT_PROC_DECL(AB2_GSRB, ab2_gsrb)
    ( const int tlo[], const int thi[],
      const int blo[], const int bhi[],
      BL_FORT_FAB_ARG(soln),
      const BL_FORT_FAB_ARG(res),
      const Real* alpha, const Real* beta,
      const BL_FORT_FAB_ARG(a),
      const BL_FORT_FAB_ARG(bx),
      const BL_FORT_FAB_ARG(by),
#if (BL_SPACEDIM == 3)
      const BL_FORT_FAB_ARG(bz),
#endif
      const BL_FORT_FAB_ARG(f0),const BL_FORT_IFAB_ARG(m0),
      const BL_FORT_FAB_ARG(f1),const BL_FORT_IFAB_ARG(m1),
      const BL_FORT_FAB_ARG(f2),const BL_FORT_IFAB_ARG(m2),
      const BL_FORT_FAB_ARG(f3),const BL_FORT_IFAB_ARG(m3),
#if (BL_SPACEDIM == 3)
      const BL_FORT_FAB_ARG(f4),const BL_FORT_IFAB_ARG(m4),
      const BL_FORT_FAB_ARG(f5),const BL_FORT_IFAB_ARG(m5),
#endif
      const int* nc, const Real dx[], const int* redBlackFlag);

BL_FORT_PROC_DECL(AB2_JACOBI, ab2_jacobi)
    ( const int lo[], const int hi[],
      BL_FORT_FAB_ARG(soln),
      const BL_FORT_FAB_ARG(res),
      const Real* alpha, const Real* beta,
      const BL_FORT_FAB_ARG(a),
      const BL_FORT_FAB_ARG(bx),
      const BL_FORT_FAB_ARG(by),
#if (BL_SPACEDIM == 3)
      const BL_FORT_FAB_ARG(bz),
#endif
      const BL_FORT_FAB_ARG(f0),const BL_FORT_IFAB_ARG(m0),
      const BL_FORT_FAB_ARG(f1),const BL_FORT_IFAB_ARG(m1),
      const BL_FORT_FAB_ARG(f2),const BL_FORT_IFAB_ARG(m2),
      const BL_FORT_FAB_ARG(f3),const BL_FORT_IFAB_ARG(m3),
#if (BL_SPACEDIM == 3)
      const BL_FORT_FAB_ARG(f4),const BL_FORT_IFAB_ARG(m4),
      const BL_FORT_FAB_ARG(f5),const BL_FORT_IFAB_ARG(m5),
#endif
      const int* nc, const Real dx[]);

#endif /*_ABec2_F_H_*/
