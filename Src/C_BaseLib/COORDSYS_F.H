#ifndef _COORDSYS_F_H_
#define _COORDSYS_F_H_

/*
** $Id: COORDSYS_F.H,v 1.5 1998-10-07 21:18:23 vince Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_SETVOL      setvol
#    define FORT_SETVOLPT      setvolpt
#    define FORT_SETAREA     setarea
#    define FORT_SETDLOGA    setdloga
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_SETVOL      SETVOL
#    define FORT_SETVOLPT      SETVOLPT
#    define FORT_SETAREA     SETAREA
#    define FORT_SETDLOGA    SETDLOGA
#  else
#  ifdef BL_AIX
#    define FORT_SETVOL      setvol
#    define FORT_SETVOLPT      setvolpt
#    define FORT_SETAREA     setarea
#    define FORT_SETDLOGA    setdloga
#  else
#    define FORT_SETVOL      setvol_
#    define FORT_SETVOLPT      setvolpt_
#    define FORT_SETAREA     setarea_
#    define FORT_SETDLOGA    setdloga_
#  endif
#  endif

#include "ArrayLim.H"

extern "C"
{
     void FORT_SETVOL (Real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       const Real* offset,
                       const Real* dx, const int* coord);

     void FORT_SETVOLPT (REAL* vol, ARLIM_P(vollo), ARLIM_P(volhi),
                        REAL* rolo, ARLIM_P(rolo), ARLIM_P(rohi),
                        REAL* rilo, ARLIM_P(rilo), ARLIM_P(rihi),
                       const REAL* dx, const int* coord);

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
}

#endif

#endif /*_COORDSYS_F_H_*/
