#ifndef _FLUXREG_F_H_
#define _FLUXREG_F_H_  "%W% %G%"

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
    void FORT_FRFINEADD (amrex::Real* reg,       ARLIM_P(rlo), ARLIM_P(rhi),
                         const amrex::Real* flx, ARLIM_P(flo), ARLIM_P(fhi),
                         const int* numcomp, const int* dir,
                         const int* ratio, const amrex::Real* mult);

    void FORT_FRFAADD (amrex::Real* reg,        ARLIM_P(rlo), ARLIM_P(rhi),
                       const amrex::Real* flx,  ARLIM_P(flo), ARLIM_P(fhi),
                       const amrex::Real* area, ARLIM_P(alo), ARLIM_P(ahi),
                       const int* numcomp, const int* dir,
                       const int* ratio, const amrex::Real* mult);

    void FORT_FRREFLUX (const int* lo, const int* hi,
			      amrex::Real* s, const int* slo, const int*shi,
			const amrex::Real* f, const int* flo, const int*fhi,
			const amrex::Real* v, const int* vlo, const int*vhi,
			const int* n, const amrex::Real* mult, 
			const int* dir, const int* islo);
#ifdef __cplusplus
}
#endif
#endif

#endif
