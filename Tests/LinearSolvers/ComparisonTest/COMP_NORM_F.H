#ifndef _COMP_NORM_F_H_
#define _COMP_NORM_F_H_

#include <AMReX_BLFort.H>

BL_FORT_PROC_DECL(LST_COMP_NORM, lst_comp_norm)
  (const int* lo, const int* hi,
   const BL_FORT_FAB_ARG(soln),
   const BL_FORT_FAB_ARG(exac),
   const BL_FORT_FAB_ARG(mask),
   const BL_FORT_FAB_ARG(volb),
   amrex_real* norm2, amrex_real* norm0,
   amrex_real* volume, int* nsoln);

#endif
