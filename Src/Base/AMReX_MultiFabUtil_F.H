#ifndef AMREX_MultiFabUtil_F_H_
#define AMREX_MultiFabUtil_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void amrex_compute_divergence (const int* lo, const int* hi,
                                   amrex_real* divu, const int* dlo, const int* dhi,
                                   const amrex_real* u, const int* ulo, const int* uhi,
#if (AMREX_SPACEDIM > 1)
                                   const amrex_real* v, const int* vlo, const int* vhi,
#endif
#if (AMREX_SPACEDIM == 3)
                                   const amrex_real* w, const int* wlo, const int* whi,
#endif
                                   const amrex_real* dxinv);
    
#ifdef __cplusplus
}
#endif


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
