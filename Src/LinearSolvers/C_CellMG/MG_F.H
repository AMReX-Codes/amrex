#ifndef _MG_F_H_
#define _MG_F_H_

/*
** $Id: MG_F.H,v 1.3 2002-04-26 22:10:02 lijewski Exp $
*/

#include <REAL.H>

#if        defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 2) 
#define FORT_AVERAGE   average2dgen
#define FORT_INTERP    interp2dgen
#endif

#if (BL_SPACEDIM == 3) 
#define FORT_AVERAGE   average3dgen
#define FORT_INTERP    interp3dgen
#endif

#else

#if (BL_SPACEDIM == 2)

#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_AVERAGE   AVERAGE2DGEN
#define FORT_INTERP    INTERP2DGEN
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_AVERAGE   average2dgen
#define FORT_INTERP    interp2dgen
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_AVERAGE   average2dgen_
#define FORT_INTERP    interp2dgen_
#endif

#endif

#if (BL_SPACEDIM == 3)

#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_AVERAGE   AVERAGE3DGEN
#define FORT_INTERP    INTERP3DGEN
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_AVERAGE   average3dgen
#define FORT_INTERP    interp3dgen
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_AVERAGE   average3dgen_
#define FORT_INTERP    interp3dgen_
#endif

#endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_AVERAGE (
        Real* crse,       ARLIM_P(crse_lo), ARLIM_P(crse_hi),
        const Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
        const int *tlo, const int *thi,
        const int *nc);

    void FORT_INTERP (
        Real* fine,       ARLIM_P(fine_lo), ARLIM_P(fine_hi),
        const Real* crse, ARLIM_P(crse_lo), ARLIM_P(crse_hi),
        const int *tlo, const int *thi,
        const int *nc);
}
#endif

#endif /*_MG_F_H_*/
