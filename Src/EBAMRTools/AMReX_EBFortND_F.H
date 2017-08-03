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
                                 const int*  refboxlo, const int* refboxhi, 
                                 const int* refrat, const int* isrc, const int* idst, const int* ncomp);

    void ebfnd_divflux(amrex::Real*       divf, const int* divflo, const int* divfhi, const int* incodivf,
                       const amrex::Real* flx0, const int* flx0lo, const int* flx0hi, const int* incoflx0,
                       const amrex::Real* flx1, const int* flx1lo, const int* flx1hi, const int* incoflx1,
                       const amrex::Real* flx2, const int* flx2lo, const int* flx2hi, const int* incoflx2,
                       const int*  gridlo, const int* gridhi, const int* multibyarea,
                       const amrex::Real* dx, const int* isrc, const int* idst, const int* ncomp);


    void ebfnd_gradlim(amrex::Real*       gph, const int* gphlo, const int* gphhi, const int* incogph,
                       const amrex::Real* phi, const int* philo, const int* phihi, const int* incophi,
                       const int*  gridlo, const int* gridhi, 
                       const int*  ncomp , const amrex::Real* dx, const int* dolimiting);

#ifdef __cplusplus
  }
#endif


  

#endif
