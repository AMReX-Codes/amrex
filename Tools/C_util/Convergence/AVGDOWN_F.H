#ifndef _AVGDOWN_F_H
#define _AVGDOWN_F_H

#ifdef BL_LANG_FORT
#    define FORT_AVGDOWN     avgdown
#    define FORT_CV_AVGDOWN  cvavgdwn
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_AVGDOWN     AVGDOWN
#    define FORT_CV_AVGDOWN  CVAVGDWN
#  else
#    define FORT_AVGDOWN     avgdown_
#    define FORT_CV_AVGDOWN  cvavgdwn_
#  endif

#include <ArrayLim.H>

extern "C" {
    
    void FORT_CV_AVGDOWN (const Real* crse, ARLIM_P(clo),     ARLIM_P(chi),
			  const int* nvar,
			  const Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
			  const int* lo, const int* hi, const int* ratio);
    
    void FORT_AVGDOWN (const Real* crse, ARLIM_P(clo),     ARLIM_P(chi),
                       const int* nvar,
                       const Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
                       const Real* cvol, ARLIM_P(cv_lo),   ARLIM_P(cv_hi), 
                       const Real* fvol, ARLIM_P(fv_lo),   ARLIM_P(fv_hi),
                       const int* lo, const int* hi, const int* ratio);
    
}
    
#endif
#endif /*_AVGDOWN_F_H_*/
