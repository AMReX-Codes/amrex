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

#include <ArrayLim.H>

extern "C"
{
    void FORT_FRFINEADD (Real* reg,       ARLIM_P(rlo), ARLIM_P(rhi),
                         const Real* flx, ARLIM_P(flo), ARLIM_P(fhi),
                         const int* numcomp, const int* dir,
                         const int* ratio, const Real* mult);

    void FORT_FRFAADD (Real* reg,        ARLIM_P(rlo), ARLIM_P(rhi),
                       const Real* flx,  ARLIM_P(flo), ARLIM_P(fhi),
                       const Real* area, ARLIM_P(alo), ARLIM_P(ahi),
                       const int* numcomp, const int* dir,
                       const int* ratio, const Real* mult);

    void FORT_FRREFLUX (Real* s,         ARLIM_P(slo), ARLIM_P(shi),
			const Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
			const Real* reg, ARLIM_P(rlo), ARLIM_P(rhi),
			const int* lo, const int* hi,
			const int* shft,
			const int* numcomp, const Real* mult);
}
#endif

#endif
