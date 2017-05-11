#ifndef _EBFriggingFortran_F_H___
#define _EBFriggingFortran_F_H___

#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>

namespace amrex
{
  extern "C"
  {
    void ebfnd_average(Real*       coar, const int* coarlo, const int* coarhi, const int* incocoar,
                       const Real* fine, const int* finelo, const int* finehi, const int* incofine,
                       const int*  coarboxlo, const int* coarboxhi, 
                       const int*   refboxlo, const int* refboxhi, 
                       const int* refrat, const int* isrc, const int* idst, const int* ncomp);


    void ebfnd_pwcinterp(Real*       fine, const int* finelo, const int* finehi, const int* incofine,
                         const Real* coar, const int* coarlo, const int* coarhi, const int* incocoar,
                         const int*  fineboxlo, const int* fineboxhi, 
                         const int* refrat, const int* isrc, const int* idst, const int* ncomp);

  }
}
#endif
