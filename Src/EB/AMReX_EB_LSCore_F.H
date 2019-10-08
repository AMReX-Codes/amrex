#ifndef AMREX_EB_LSCore_F_H_
#define AMREX_EB_LSCore_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    // these are interfaces to fortran subroutines

    void amrex_eb_levelset_error (int * tag, const int * tag_lo, const int * tag_hi,
                                  const BL_FORT_FAB_ARG_3D(state),
                                  const int * tagval, const int * clearval,
                                  const int * lo,     const int * hi,
                                  const amrex_real * dx,   const amrex_real * problo,
                                  const amrex_real * time, const amrex_real * phierr);

    void amrex_eb_volfrac_error (int * tag, const int * tag_lo, const int * tag_hi,
                                 const BL_FORT_FAB_ARG_3D(state),
                                 const int * tagval, const int * clearval, const double * tol,
                                 const int * lo,     const int * hi,
                                 const amrex_real * dx,   const amrex_real * problo);

    void amrex_eb_phifill (BL_FORT_FAB_ARG_3D(state),
                           const int* dlo, const int* dhi,
                           const amrex_real* dx, const amrex_real* glo,
                           const amrex_real* time, const int* bc);

#ifdef __cplusplus
}
#endif

#endif
