#ifndef _COORDSYS_F_H_
#define _COORDSYS_F_H_  

#ifdef BL_LANG_FORT
#    define FORT_SETVOL      setvol
#    define FORT_SETAREA     setarea
#    define FORT_SETDLOGA    setdloga
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_SETVOL      SETVOL
#    define FORT_SETAREA     SETAREA
#    define FORT_SETDLOGA    SETDLOGA
#  else
#    define FORT_SETVOL      setvol_
#    define FORT_SETAREA     setarea_
#    define FORT_SETDLOGA    setdloga_
#  endif

#include "ArrayLim.H"

  extern "C" {
     void FORT_SETVOL (Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       const Real* offset,
                       const Real* dx, const int* coord);

#if (BL_SPACEDIM == 2)
     void FORT_SETDLOGA (Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
			 const Real* offset,
			 const Real* dx, const int* dir,
			 const int* coord);
#endif

     void FORT_SETAREA (Real* area,  ARLIM_P(alo), ARLIM_P(ahi),
                        const Real* offset,
                        const Real* dx, const int* dir,
                        const int* coord);
  };

#endif

#endif
