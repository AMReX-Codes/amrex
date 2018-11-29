#ifndef AMREX_MultiFabUtil_F_H_
#define AMREX_MultiFabUtil_F_H_

#include <AMReX_BLFort.H>


BL_FORT_PROC_DECL(BL_AVGDOWN_WITH_VOL, bl_avgdown_with_vol)
  (const int* ovlp_lo, const int* ovlp_hi,
   const BL_FORT_FAB_ARG(fine),
   BL_FORT_FAB_ARG(crse),
   const BL_FORT_FAB_ARG(fvolume),
   const int* ratio, const int* ncomp);

#endif
