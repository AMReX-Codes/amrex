#ifndef _INTERPBNDRY_F_H_
#define _INTERPBNDRY_F_H_

/*
** $Id: INTERPBNDRYDATA_F.H,v 1.1 1998-03-24 07:06:07 almgren Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_BDINTERPXLO  bdintrpxlo
#    define FORT_BDINTERPXHI  bdintrpxhi
#  if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  bdintrpylo
#    define FORT_BDINTERPYHI  bdintrpyhi
#  endif
#  if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  bdintrpzlo
#    define FORT_BDINTERPZHI  bdintrpzhi
#  endif
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_BDINTERPXLO  BDINTRPXLO
#    define FORT_BDINTERPXHI  BDINTRPXHI
#    if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  BDINTRPYLO
#    define FORT_BDINTERPYHI  BDINTRPYHI
#    endif
#    if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  BDINTRPZLO
#    define FORT_BDINTERPZHI  BDINTRPZHI
#    endif
#  else
#    define FORT_BDINTERPXLO  bdintrpxlo_
#    define FORT_BDINTERPXHI  bdintrpxhi_
#    if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  bdintrpylo_
#    define FORT_BDINTERPYHI  bdintrpyhi_
#    endif
#    if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  bdintrpzlo_
#    define FORT_BDINTERPZHI  bdintrpzhi_
#    endif
#endif

#include <ArrayLim.H>

typedef void BDInterpFunc(Real* bndr, ARLIM_P(blo), ARLIM_P(bhi),
                          const int* lo, const int* hi,
                          ARLIM_P(cblo), ARLIM_P(cbhi),
                          const int* nvar, const int* ratio,
                          const int* not_covered,
                          const int* mask,  ARLIM_P(mlo), ARLIM_P(mhi),
                          const Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                          Real* derives);
extern "C"
{
    BDInterpFunc FORT_BDINTERPXLO;
    BDInterpFunc FORT_BDINTERPXHI;

#if (BL_SPACEDIM > 1)    
    BDInterpFunc FORT_BDINTERPYLO;
    BDInterpFunc FORT_BDINTERPYHI;
#endif    

#if (BL_SPACEDIM > 2)    
    BDInterpFunc FORT_BDINTERPZLO;
    BDInterpFunc FORT_BDINTERPZHI;
#endif    

}
#endif

#endif /*_INTERPBNDRY_F_H_*/
