#ifndef _FLUXREG_F_H_
#define _FLUXREG_F_H_

/*
** $Id: FLUXREG_F.H,v 1.4 1997-12-14 23:34:28 lijewski Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_FRCRSEINIT  frcrseinit
#    define FORT_FRFINEADD   frfineadd
#    define FORT_FRCAINIT    frcainit
#    define FORT_FRFAADD     frfaadd
#    define FORT_FRREFLUX    frreflux
#    define FORT_FRCVREFLUX  frcvreflux
#    define FORT_SCALADDTO   scaladdto
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_FRCRSEINIT  FRCRSEINIT
#    define FORT_FRFINEADD   FRFINEADD
#    define FORT_FRCAINIT    FRCAINIT
#    define FORT_FRFAADD     FRFAADD
#    define FORT_FRREFLUX    FRREFLUX
#    define FORT_FRCVREFLUX  FRCVREFLUX
#    define FORT_SCALADDTO   SCALADDTO
#  else
#    define FORT_FRCRSEINIT  frcrseinit_
#    define FORT_FRFINEADD   frfineadd_
#    define FORT_FRCAINIT    frcainit_
#    define FORT_FRFAADD     frfaadd_
#    define FORT_FRREFLUX    frreflux_
#    define FORT_FRCVREFLUX  frcvreflux_
#    define FORT_SCALADDTO   scaladdto_
#  endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_FRCRSEINIT (Real* reg,       ARLIM_P(rlo), ARLIM_P(rhi),
                          const Real* flx, ARLIM_P(flo), ARLIM_P(fhi),
                          const int* lo, const int* hi,
                          const int* numcomp, const int* dir,
                          const Real* mult);

    void FORT_FRCAINIT (Real* reg,        ARLIM_P(rlo), ARLIM_P(rhi),
                        const Real* flx,  ARLIM_P(flo), ARLIM_P(fhi),
                        const Real* area, ARLIM_P(alo), ARLIM_P(ahi),
                        const int* lo, const int* hi,
                        const int* numcomp, const int* dir,
                        const Real* mult);

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
                        const int* numcomp, const Real* mult);

    void FORT_FRCVREFLUX (Real* s, ARLIM_P(slo), ARLIM_P(shi),
                          const Real* dx,
                          const Real* reg, ARLIM_P(rlo), ARLIM_P(rhi),
                          const int* lo, const int* hi,
                          const int* numcomp, const Real* mult);

    void FORT_SCALADDTO (Real * flx,         ARLIM_P(flo), ARLIM_P(fhi),
                         const Real * area,  ARLIM_P(alo), ARLIM_P(ahi),
                         Real * reg,         ARLIM_P(rlo), ARLIM_P(rhi),
                         const int * lo, const int * hi,
                         const int * numcomp, const Real * mult);

};
#endif

#endif /*_FLUXREG_F_H_*/
