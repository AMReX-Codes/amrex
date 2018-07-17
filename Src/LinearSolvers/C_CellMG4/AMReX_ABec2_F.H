#ifndef AMREX_ABec2_F_H_
#define AMREX_ABec2_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void amrex_ab2_bndryrlx
    ( const int lo[], const int hi[],
      BL_FORT_FAB_ARG(f),
      const BL_FORT_IFAB_ARG(m),
      const int* cdir, const int* bct,
      const amrex_real* bcl, const int* maxorder, const amrex_real dx[] );

    void amrex_ab2_gsrb
    ( const int tlo[], const int thi[],
      const int blo[], const int bhi[],
      BL_FORT_FAB_ARG(soln),
      const BL_FORT_FAB_ARG(res),
      const amrex_real* alpha, const amrex_real* beta,
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
      const int* nc, const amrex_real dx[], const int* redBlackFlag);

    void amrex_ab2_jacobi
    ( const int lo[], const int hi[],
      BL_FORT_FAB_ARG(soln),
      const BL_FORT_FAB_ARG(res),
      const amrex_real* alpha, const amrex_real* beta,
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
      const int* nc, const amrex_real dx[]);

#ifdef __cplusplus
}
#endif

#endif /*_ABec2_F_H_*/
