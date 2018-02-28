
#ifndef _INTERPBNDRY_F_H_
#define _INTERPBNDRY_F_H_

#ifdef BL_LANG_FORT
#    define FORT_BDINTERPXLO  bdintrpxlo
#    define FORT_BDINTERPXHI  bdintrpxhi
#  if (AMREX_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  bdintrpylo
#    define FORT_BDINTERPYHI  bdintrpyhi
#  endif
#  if (AMREX_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  bdintrpzlo
#    define FORT_BDINTERPZHI  bdintrpzhi
#  endif
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_BDINTERPXLO  BDINTRPXLO
#    define FORT_BDINTERPXHI  BDINTRPXHI
#    if (AMREX_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  BDINTRPYLO
#    define FORT_BDINTERPYHI  BDINTRPYHI
#    endif
#    if (AMREX_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  BDINTRPZLO
#    define FORT_BDINTERPZHI  BDINTRPZHI
#    endif
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_BDINTERPXLO  bdintrpxlo
#    define FORT_BDINTERPXHI  bdintrpxhi
#    if (AMREX_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  bdintrpylo
#    define FORT_BDINTERPYHI  bdintrpyhi
#    endif
#    if (AMREX_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  bdintrpzlo
#    define FORT_BDINTERPZHI  bdintrpzhi
#    endif
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_BDINTERPXLO  bdintrpxlo_
#    define FORT_BDINTERPXHI  bdintrpxhi_
#    if (AMREX_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  bdintrpylo_
#    define FORT_BDINTERPYHI  bdintrpyhi_
#    endif
#    if (AMREX_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  bdintrpzlo_
#    define FORT_BDINTERPZHI  bdintrpzhi_
#    endif
#endif

#include <AMReX_ArrayLim.H>
#include <AMReX_REAL.H>

typedef void BDInterpFunc(amrex_real* bndr, ARLIM_P(blo), ARLIM_P(bhi),
                          const int* lo, const int* hi,
                          ARLIM_P(cblo), ARLIM_P(cbhi),
                          const int* nvar, const int* ratio,
                          const int* not_covered,
                          const int* mask,  ARLIM_P(mlo), ARLIM_P(mhi),
                          const amrex_real* crse, ARLIM_P(clo), ARLIM_P(chi),
                          amrex_real* derives, int* max_order);

#ifdef __cplusplus
extern "C"
{
#endif
    BDInterpFunc FORT_BDINTERPXLO;
    BDInterpFunc FORT_BDINTERPXHI;

#if (AMREX_SPACEDIM > 1)    
    BDInterpFunc FORT_BDINTERPYLO;
    BDInterpFunc FORT_BDINTERPYHI;
#endif    

#if (AMREX_SPACEDIM > 2)    
    BDInterpFunc FORT_BDINTERPZLO;
    BDInterpFunc FORT_BDINTERPZHI;
#endif    

#ifdef __cplusplus
}
#endif
#endif

#endif /*_INTERPBNDRY_F_H_*/
