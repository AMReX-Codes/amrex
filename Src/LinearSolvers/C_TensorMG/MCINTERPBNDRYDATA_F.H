#ifndef _INTERPBNDRY_F_H_
#define _INTERPBNDRY_F_H_  "%W% %G%"

#ifdef BL_LANG_FORT
#    define FORT_BDINTERPXLO  mcbdintrpxlo
#    define FORT_BDINTERPXHI  mcbdintrpxhi
#    define FORT_BDIDERIVXLO  mcbdiderivxlo
#    define FORT_BDIDERIVXHI  mcbdiderivxhi
#  if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  mcbdintrpylo
#    define FORT_BDINTERPYHI  mcbdintrpyhi
#    define FORT_BDIDERIVYLO  mcbdiderivylo
#    define FORT_BDIDERIVYHI  mcbdiderivyhi
#  endif
#  if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  mcbdintrpzlo
#    define FORT_BDINTERPZHI  mcbdintrpzhi
#    define FORT_BDIDERIVZLO  mcbdiderivzlo
#    define FORT_BDIDERIVZHI  mcbdiderivzhi
#  endif
#endif

#ifdef BL_LANG_CC
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_BDINTERPXLO  MCBDINTRPXLO
#    define FORT_BDINTERPXHI  MCBDINTRPXHI
#    define FORT_BDIDERIVXLO  MCBDIDERIVXLO
#    define FORT_BDIDERIVXHI  MCBDIDERIVXHI
#    if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  MCBDINTRPYLO
#    define FORT_BDINTERPYHI  MCBDINTRPYHI
#    define FORT_BDIDERIVYLO  MCBDIDERIVYLO
#    define FORT_BDIDERIVYHI  MCBDIDERIVYHI
#    endif
#    if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  MCBDINTRPZLO
#    define FORT_BDINTERPZHI  MCBDINTRPZHI
#    define FORT_BDIDERIVZLO  MCBDIDERIVZLO
#    define FORT_BDIDERIVZHI  MCBDIDERIVZHI
#    endif
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_BDINTERPXLO  mcbdintrpxlo
#    define FORT_BDINTERPXHI  mcbdintrpxhi
#    define FORT_BDIDERIVXLO  mcbdiderivxlo
#    define FORT_BDIDERIVXHI  mcbdiderivxhi
#    if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  mcbdintrpylo
#    define FORT_BDINTERPYHI  mcbdintrpyhi
#    define FORT_BDIDERIVYLO  mcbdiderivylo
#    define FORT_BDIDERIVYHI  mcbdiderivyhi
#    endif
#    if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  mcbdintrpzlo
#    define FORT_BDINTERPZHI  mcbdintrpzhi
#    define FORT_BDIDERIVZLO  mcbdiderivzlo
#    define FORT_BDIDERIVZHI  mcbdiderivzhi
#    endif
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_BDINTERPXLO  mcbdintrpxlo_
#    define FORT_BDINTERPXHI  mcbdintrpxhi_
#    define FORT_BDIDERIVXLO  mcbdiderivxlo_
#    define FORT_BDIDERIVXHI  mcbdiderivxhi_
#    if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  mcbdintrpylo_
#    define FORT_BDINTERPYHI  mcbdintrpyhi_
#    define FORT_BDIDERIVYLO  mcbdiderivylo_
#    define FORT_BDIDERIVYHI  mcbdiderivyhi_
#    endif
#    if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  mcbdintrpzlo_
#    define FORT_BDINTERPZHI  mcbdintrpzhi_
#    define FORT_BDIDERIVZLO  mcbdiderivzlo_
#    define FORT_BDIDERIVZHI  mcbdiderivzhi_
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
			  Real* derives, const Real* hfine);
typedef void BDPhysDerivative(Real* bndr, ARLIM_P(blo), ARLIM_P(bhi),
			      const int* lo, const int* hi,
			      const Real* finefab, ARLIM_P(flo), ARLIM_P(fhi),
			      const int* nvar, const Real* hfine);
extern "C"
{
    BDInterpFunc FORT_BDINTERPXLO;
    BDInterpFunc FORT_BDINTERPXHI;
    BDPhysDerivative FORT_BDIDERIVXLO;
    BDPhysDerivative FORT_BDIDERIVXHI;

#if (BL_SPACEDIM > 1)    
    BDInterpFunc FORT_BDINTERPYLO;
    BDInterpFunc FORT_BDINTERPYHI;
    BDPhysDerivative FORT_BDIDERIVYLO;
    BDPhysDerivative FORT_BDIDERIVYHI;
#endif    

#if (BL_SPACEDIM > 2)    
    BDInterpFunc FORT_BDINTERPZLO;
    BDInterpFunc FORT_BDINTERPZHI;
    BDPhysDerivative FORT_BDIDERIVZLO;
    BDPhysDerivative FORT_BDIDERIVZHI;
#endif    

}
#endif

#endif
