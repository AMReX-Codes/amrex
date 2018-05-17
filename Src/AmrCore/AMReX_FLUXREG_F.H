#ifndef AMREX_FLUXREG_F_H_
#define AMREX_FLUXREG_F_H_

#ifdef BL_LANG_FORT
#    define FORT_FRFINEADD   frfineadd
#    define FORT_FRFAADD     frfaadd
#    define FORT_FRREFLUX    frreflux
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_FRFINEADD   FRFINEADD
#    define FORT_FRFAADD     FRFAADD
#    define FORT_FRREFLUX    FRREFLUX
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_FRFINEADD   frfineadd
#    define FORT_FRFAADD     frfaadd
#    define FORT_FRREFLUX    frreflux
#  else
#    define FORT_FRFINEADD   frfineadd_
#    define FORT_FRFAADD     frfaadd_
#    define FORT_FRREFLUX    frreflux_
#  endif
#  endif

#include <AMReX_ArrayLim.H>

#include <AMReX_REAL.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void FORT_FRFINEADD (amrex_real* reg,       ARLIM_P(rlo), ARLIM_P(rhi),
                         const amrex_real* flx, ARLIM_P(flo), ARLIM_P(fhi),
                         const int* numcomp, const int* dir,
                         const int* ratio, const amrex_real* mult);

    void FORT_FRFAADD (amrex_real* reg,        ARLIM_P(rlo), ARLIM_P(rhi),
                       const amrex_real* flx,  ARLIM_P(flo), ARLIM_P(fhi),
                       const amrex_real* area, ARLIM_P(alo), ARLIM_P(ahi),
                       const int* numcomp, const int* dir,
                       const int* ratio, const amrex_real* mult);

    void FORT_FRREFLUX (const int* lo, const int* hi,
			      amrex_real* s, const int* slo, const int*shi,
			const amrex_real* f, const int* flo, const int*fhi,
			const amrex_real* v, const int* vlo, const int*vhi,
			const int* n, const amrex_real* mult, 
			const int* dir, const int* islo);
#ifdef __cplusplus
}
#endif
#endif

#endif
