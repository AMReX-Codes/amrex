#ifndef _INTERP_NAMES_H
#define _INTERP_NAMES_H

#ifdef BL_LANG_FORT
#    define FORT_NBINTERP    nbinterp
#    define FORT_CBINTERP    cbinterp
#    define FORT_CCINTERP    ccinterp
#    define FORT_CCINTERP2   ccinterp2
#    define FORT_PCINTERP    pcinterp
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_NBINTERP    NBINTERP
#    define FORT_CBINTERP    CBINTERP
#    define FORT_CCINTERP    CCINTERP
#    define FORT_CCINTERP2   CCINTERP2
#    define FORT_PCINTERP    PCINTERP
#  else
#    define FORT_NBINTERP    nbinterp_
#    define FORT_CBINTERP    cbinterp_
#    define FORT_CCINTERP    ccinterp_
#    define FORT_CCINTERP2   ccinterp2_
#    define FORT_PCINTERP    pcinterp_
#  endif

#include <ArrayLim.H>

  extern "C" {
     void FORT_NBINTERP (const REAL* crse, ARLIM_P(clo), ARLIM_P(chi),
                         ARLIM_P(cblo), ARLIM_P(cbhi),
                         REAL* fine, ARLIM_P(flo), ARLIM_P(fhi),
                         ARLIM_P(fblo), ARLIM_P(fbhi),
			 const int* lratio, const int* nvar,
			 REAL* slope, const int* num_slope,
			 REAL* strip, const int* strip_lo, const int* strip_hi);
      
     void FORT_CBINTERP (const REAL* crse, ARLIM_P(clo), ARLIM_P(chi),
                         ARLIM_P(cblo), ARLIM_P(cbhi),
                         REAL* fine, ARLIM_P(flo), ARLIM_P(fhi),
                         ARLIM_P(fblo), ARLIM_P(fbhi),
			 const int* lratio, const int* nvar,
			 REAL* slope, const int* num_slope,
			 REAL* strip, const int* strip_lo, const int* strip_hi);

     void FORT_CCINTERP (REAL* fine, ARLIM_P(flo), ARLIM_P(fhi),
                         ARLIM_P(fblo), ARLIM_P(fbhi),
			 const int* nvar, const int* lratio,
			 const REAL* crse, const int* clo, const int* chi,
                         ARLIM_P(cblo), ARLIM_P(cbhi),
     		         const int* fslo, const int* fshi,
			 REAL* cslope, const int* clen,
			 REAL* fslope, REAL* fdat, const int* flen,
			 REAL* voff,
			 const int* bc, const int* lim_slope,
		 D_DECL(const REAL* fvcx,const REAL* fvcy, const REAL* fvcz),
		 D_DECL(const REAL* cvcx,const REAL* cvcy, const REAL* cvcz));

     void FORT_CCINTERP2 (REAL* fine, const int* flo, const int* fhi,
                          ARLIM_P(fblo), ARLIM_P(fbhi),
			  const int* nvar, const int* lratio,
			  const REAL* crse, const int* clo, const int* chi,
                          ARLIM_P(cblo), ARLIM_P(cbhi),
     		          const int* fslo, const int* fshi,
			  REAL* cslope, const int* clen,
			  REAL* fslope, REAL* fdat, const int* flen,
			  REAL* voff,
			  const int* bc, const int* lim_slope,
		  D_DECL(const REAL* fvcx,const REAL* fvcy, const REAL* fvcz),
		  D_DECL(const REAL* cvcx,const REAL* cvcy, const REAL* cvcz));

     void FORT_PCINTERP (const REAL* crse, ARLIM_P(clo), ARLIM_P(chi),
                         const int* cblo, const int* cbhi,
			 REAL*  fine, ARLIM_P(flo), ARLIM_P(fhi),
			 const int* fblo, const int* fbhi,
			 const int* longdir, const int* ratio,
			 const int* nvar, REAL* strip,
			 const int* strip_lo, const int* strip_hi);
  };
#endif

#endif
