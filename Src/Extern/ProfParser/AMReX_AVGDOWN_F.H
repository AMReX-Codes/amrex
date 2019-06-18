#ifndef AMREX_AVGDOWN_F_H
#define AMREX_AVGDOWN_F_H

#ifdef BL_LANG_FORT
#    define FORT_AVGDOWN         avgdown
#    define FORT_MAXVAL_AVGDOWN  maxvalavgdwn
#    define FORT_CV_AVGDOWN      cvavgdwn
#    define FORT_CV_AVGDOWN_STAG cvavgdwnstag
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_AVGDOWN         AVGDOWN
#    define FORT_MAXVAL_AVGDOWN  MAXVALAVGDWN
#    define FORT_CV_AVGDOWN      CVAVGDWN
#    define FORT_CV_AVGDOWN_STAG CVAVGDWNSTAG
#  else
#    define FORT_AVGDOWN         avgdown_
#    define FORT_MAXVAL_AVGDOWN  maxvalavgdwn_
#    define FORT_CV_AVGDOWN      cvavgdwn_
#    define FORT_CV_AVGDOWN_STAG cvavgdwnstag_
#  endif

#include <AMReX_ArrayLim.H>

using amrex::Real;

extern "C" {
    
    void FORT_MAXVAL_AVGDOWN (const Real* crse, AMREX_ARLIM_P(clo),     AMREX_ARLIM_P(chi),
			  const int* nvar,
			  const Real* fine, AMREX_ARLIM_P(fine_lo), AMREX_ARLIM_P(fine_hi),
			  const int* lo, const int* hi, const int* ratio);

    void FORT_CV_AVGDOWN (const Real* crse, AMREX_ARLIM_P(clo),     AMREX_ARLIM_P(chi),
			  const int* nvar,
			  const Real* fine, AMREX_ARLIM_P(fine_lo), AMREX_ARLIM_P(fine_hi),
			  const int* lo, const int* hi, const int* ratio);

    void FORT_CV_AVGDOWN_STAG (const int* nodal_dir, 
  			       const Real* crse, AMREX_ARLIM_P(clo),     AMREX_ARLIM_P(chi),
    			       const int* nvar,
			       const Real* fine, AMREX_ARLIM_P(fine_lo), AMREX_ARLIM_P(fine_hi),
			       const int* lo, const int* hi, const int* ratio);
    
    void FORT_AVGDOWN (const Real* crse, AMREX_ARLIM_P(clo),     AMREX_ARLIM_P(chi),
                       const int* nvar,
                       const Real* fine, AMREX_ARLIM_P(fine_lo), AMREX_ARLIM_P(fine_hi),
                       const Real* cvol, AMREX_ARLIM_P(cv_lo),   AMREX_ARLIM_P(cv_hi), 
                       const Real* fvol, AMREX_ARLIM_P(fv_lo),   AMREX_ARLIM_P(fv_hi),
                       const int* lo, const int* hi, const int* ratio);
    
}
    
#endif
#endif /*_AVGDOWN_F_H_*/
