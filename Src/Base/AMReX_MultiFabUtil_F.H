#ifndef _MtltiFabUtil_F_H_
#define _MtltiFabUtil_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void amrex_fort_avg_nd_to_cc(const int* lo, const int* hi, const int* ncomp,
                                 BL_FORT_FAB_ARG(cc), const BL_FORT_FAB_ARG(nd));

    void amrex_fort_int_to_real (const int* lo, const int* hi, const int* ncomp,
                                 amrex_real* rdata, const int* rlo, const int* rhi,
                                 const int* idata, const int* ilo, const int* ihi);

    void amrex_fill_slice(const int* lo, const int* hi,
                          const amrex_real* full_data, const int* flo, const int* fhi,
                          amrex_real* slice_data, const int* slo, const int* shi,
                          const int* fstart, const int* nfull, const int* ncomp);

#ifdef __cplusplus
}
#endif

BL_FORT_PROC_DECL(BL_AVG_EG_TO_CC, bl_avg_eg_to_cc)
  (const int* lo, const int* hi,
   BL_FORT_FAB_ARG(cc),
   AMREX_D_DECL(const BL_FORT_FAB_ARG(Ex),
	  const BL_FORT_FAB_ARG(Ey),
	  const BL_FORT_FAB_ARG(Ez)));

BL_FORT_PROC_DECL(BL_AVG_FC_TO_CC, bl_avg_fc_to_cc)
  (const int* lo, const int* hi,
   BL_FORT_FAB_ARG(cc),
   AMREX_D_DECL(const BL_FORT_FAB_ARG(fx),
	  const BL_FORT_FAB_ARG(fy),
	  const BL_FORT_FAB_ARG(fz)),
   const amrex_real* dx, const amrex_real* problo, const int& coord_type);

BL_FORT_PROC_DECL(BL_AVG_CC_TO_FC, bl_avg_cc_to_fc)
  (const int* xlo, const int* xhi,
#if (AMREX_SPACEDIM > 1)
   const int* ylo, const int* yhi,
#endif
#if (AMREX_SPACEDIM == 3)
   const int* zlo, const int* zhi,
#endif
   AMREX_D_DECL(BL_FORT_FAB_ARG(fx),
	  BL_FORT_FAB_ARG(fy),
	  BL_FORT_FAB_ARG(fz)),
   const BL_FORT_FAB_ARG(cc),
   const amrex_real* dx, const amrex_real* problo, const int& coord_type);

BL_FORT_PROC_DECL(BL_AVGDOWN_FACES, bl_avgdown_faces)
  (const int* lo, const int* hi,
   const BL_FORT_FAB_ARG(fine),
   BL_FORT_FAB_ARG(crse),
   const int* ratio, const int& n, const int& ncomp); 

BL_FORT_PROC_DECL(BL_AVGDOWN_EDGES, bl_avgdown_edges)
  (const int* lo, const int* hi,
   const BL_FORT_FAB_ARG(fine),
   BL_FORT_FAB_ARG(crse),
   const int* ratio, const int& n, const int& ncomp); 

BL_FORT_PROC_DECL(BL_AVGDOWN, bl_avgdown)
  (const int* ovlp_lo, const int* ovlp_hi,
   const BL_FORT_FAB_ARG(fine),
   BL_FORT_FAB_ARG(crse),
   const int* ratio, const int* ncomp);

BL_FORT_PROC_DECL(BL_AVGDOWN_NODES, bl_avgdown_nodes)
  (const int* ovlp_lo, const int* ovlp_hi,
   const BL_FORT_FAB_ARG(fine),
   BL_FORT_FAB_ARG(crse),
   const int* ratio, const int* ncomp);

BL_FORT_PROC_DECL(BL_AVGDOWN_WITH_VOL, bl_avgdown_with_vol)
  (const int* ovlp_lo, const int* ovlp_hi,
   const BL_FORT_FAB_ARG(fine),
   BL_FORT_FAB_ARG(crse),
   const BL_FORT_FAB_ARG(fvolume),
   const int* ratio, const int* ncomp);

#endif
