#ifndef _MtltiFabUtil_F_H_
#define _MtltiFabUtil_F_H_

#include <BLFort.H>

BL_FORT_PROC_DECL(BL_AVG_FC_TO_CC, bl_avg_fc_to_cc)
  (const int* lo, const int* hi,
   BL_FORT_FAB_ARG(cc),
   D_DECL(const BL_FORT_FAB_ARG(fx),
	  const BL_FORT_FAB_ARG(fy),
	  const BL_FORT_FAB_ARG(fz)),
   const Real* dx, const Real* problo, const int& coord_type);

BL_FORT_PROC_DECL(BL_AVGDOWN_FACES, bl_avgdown_faces)
  (const int* lo, const int* hi,
   BL_FORT_FAB_ARG(fine),
   BL_FORT_FAB_ARG(crse),
   const int* ratio, const int& n); 


#endif
