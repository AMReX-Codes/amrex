#ifndef _MFTOWER_F_H
#define _MFTOWER_F_H

#ifdef BL_LANG_FORT
#    define FORT_CC_TO_EC_GRAD    cctoecgrad
#    define FORT_EC_TO_CC_DIV     ectoccdiv
#    define FORT_COARSEN_EC       coarsenec
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_CC_TO_EC_GRAD    CCTOECGRAD
#    define FORT_EC_TO_CC_DIV     ECTOCCDIV
#    define FORT_COARSEN_EC       COARSENEC
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_CC_TO_EC_GRAD    cctoecgrad
#    define FORT_EC_TO_CC_DIV     ectoccdiv
#    define FORT_COARSEN_EC       coarsenec
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_CC_TO_EC_GRAD    cctoecgrad_
#    define FORT_EC_TO_CC_DIV     ectoccdiv_
#    define FORT_COARSEN_EC       coarsenec_
#  endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_CC_TO_EC_GRAD (Real* ecg,      ARLIM_P(elo), ARLIM_P(ehi),
                             const Real* cc, ARLIM_P(clo), ARLIM_P(chi),
                             const int*  ccvlo, const int*  ccvhi,
                             const Real* dx, const Real* mult, const int* dir,const int* nc);

    void FORT_EC_TO_CC_DIV (Real* ccg,      ARLIM_P(clo), ARLIM_P(chi),
                            const Real* ec, ARLIM_P(elo), ARLIM_P(ehi),
                            const Real* a,  ARLIM_P(alo), ARLIM_P(ahi),
                            const Real* v,  ARLIM_P(vlo), ARLIM_P(vhi),
                            const int*  ccvlo, const int*  ccvhi,
                            const Real* mult,
                            const int* dir, const int* nc);

    void FORT_COARSEN_EC(Real* cefab,       ARLIM_P(celo), ARLIM_P(cehi),
                         const Real* fefab, ARLIM_P(felo), ARLIM_P(fehi),
                         const int* elo, const int* ehi, const int* ratio, 
                         const int* dir, const int* nc);
}

#endif
#endif /*_MFTOWER_F_H*/
