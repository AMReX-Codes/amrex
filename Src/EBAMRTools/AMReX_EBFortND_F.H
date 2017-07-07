#ifndef _EBFriggingFortran_F_H___
#define _EBFriggingFortran_F_H___

#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
  extern "C"
  {
#endif
    void ebfnd_average(amrex::Real*       coar, const int* coarlo, const int* coarhi, const int* incocoar,
                       const amrex::Real* fine, const int* finelo, const int* finehi, const int* incofine,
                       const int*  coarboxlo, const int* coarboxhi, 
                       const int*   refboxlo, const int* refboxhi,
                       const int* refrat, const int* isrc, const int* idst, const int* ncomp);


    void ebfnd_average_face(amrex::Real*       coar, const int* coarlo, const int* coarhi, const int* incocoar,
                            const amrex::Real* fine, const int* finelo, const int* finehi, const int* incofine,
                            const int*  coarboxlo, const int* coarboxhi,
                            const int*  refboxlo, const int* refboxhi,
                            const int* facedir   , const int* refrat   , 
                            const int* isrc, const int* idst, const int* ncomp);


    void ebfnd_pwcinterp(amrex::Real*       fine, const int* finelo, const int* finehi, const int* incofine,
                         const amrex::Real* coar, const int* coarlo, const int* coarhi, const int* incocoar,
                         const int*  fineboxlo, const int* fineboxhi, 
                         const int* refrat, const int* isrc, const int* idst, const int* ncomp);


    void ebfnd_pwlinterp_nobound(amrex::Real*       fine, const int* finelo, const int* finehi, const int* incofine,
                                 const amrex::Real* coar, const int* coarlo, const int* coarhi, const int* incocoar,
                                 const int*  fineboxlo, const int* fineboxhi, 
                                 const int* refrat, const int* isrc, const int* idst, const int* ncomp);

    void ebfnd_pwqinterp_nobound(amrex::Real*       fine, const int* finelo, const int* finehi, const int* incofine,
                                 const amrex::Real* coar, const int* coarlo, const int* coarhi, const int* incocoar,
                                 const int*  fineboxlo, const int* fineboxhi, 
                                 const int* refrat, const int* isrc, const int* idst, const int* ncomp);


    void ebfnd_pwl_incr_at_bound(amrex::Real*       fine, const int* finelo, const int* finehi, const int* incofine,
                                 const amrex::Real* coar, const int* coarlo, const int* coarhi, const int* incocoar,
                                 const int* has_lo, const int* has_hi, const int* idir,
                                 const int*  coar_lo_box_lo, const int* coar_lo_box_hi,
                                 const int*  coar_hi_box_lo, const int* coar_hi_box_hi,
                                 const int*  coar_cen_box_lo, const int* coar_cen_box_hi,
                                 const int*   refboxlo, const int* refboxhi, 
                                 const int* refrat, const int* isrc, const int* idst, const int* ncomp);

#ifdef __cplusplus
  }
#endif


  

#endif
