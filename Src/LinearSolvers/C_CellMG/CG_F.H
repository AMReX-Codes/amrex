#ifndef _CG_F_H_
#define _CG_F_H_

/*
** $Id: CG_F.H,v 1.4 2002-04-26 22:10:02 lijewski Exp $
*/

#include <REAL.H>

#if        defined(BL_LANG_FORT)

#if (BL_SPACEDIM == 2)
#define FORT_CGUPDATE  cgupdate2dgen
#define FORT_CGADVCP   cgadvcp2dgen
#define FORT_CGXDOTY   cgxdoty2dgen
#define FORT_CGSXAY   cgsxay2dgen
#endif

#if (BL_SPACEDIM == 3)
#define FORT_CGUPDATE  cgupdate3dgen
#define FORT_CGADVCP   cgadvcp3dgen
#define FORT_CGXDOTY   cgxdoty3dgen
#define FORT_CGSXAY   cgsxay3dgen
#endif

#else

#if defined(BL_FORT_USE_UPPERCASE)
#  if (BL_SPACEDIM == 2)
#    define FORT_CGUPDATE  CGUPDATE2DGEN
#    define FORT_CGADVCP   CGADVCP2DGEN
#    define FORT_CGXDOTY   CGXDOTY2DGEN
#    define FORT_CGSXAY   CGSXAY2DGEN
#  endif
#  if (BL_SPACEDIM == 3)
#    define FORT_CGUPDATE  CGUPDATE3DGEN
#    define FORT_CGADVCP   CGADVCP3DGEN
#    define FORT_CGXDOTY   CGXDOTY3DGEN
#    define FORT_CGSXAY   CGSXAY3DGEN
#  endif
#elif defined(BL_FORT_USE_LOWERCASE)
#  if (BL_SPACEDIM == 2)
#    define FORT_CGUPDATE  cgupdate2dgen
#    define FORT_CGADVCP   cgadvcp2dgen
#    define FORT_CGXDOTY   cgxdoty2dgen
#    define FORT_CGSXAY   cgsxay2dgen
#  endif
#  if (BL_SPACEDIM == 3)
#    define FORT_CGUPDATE  cgupdate3dgen
#    define FORT_CGADVCP   cgadvcp3dgen
#    define FORT_CGXDOTY   cgxdoty3dgen
#    define FORT_CGSXAY   cgsxay3dgen
#  endif
#elif defined(BL_FORT_USE_UNDERSCORE)
#  if (BL_SPACEDIM == 2)
#    define FORT_CGUPDATE  cgupdate2dgen_
#    define FORT_CGADVCP   cgadvcp2dgen_
#    define FORT_CGXDOTY   cgxdoty2dgen_
#    define FORT_CGSXAY   cgsxay2dgen_
#  endif
#  if (BL_SPACEDIM == 3)
#    define FORT_CGUPDATE  cgupdate3dgen_
#    define FORT_CGADVCP   cgadvcp3dgen_
#    define FORT_CGXDOTY   cgxdoty3dgen_
#    define FORT_CGSXAY   cgsxay3dgen_
#  endif
#endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_CGUPDATE (
        Real* phi,  ARLIM_P(phi_lo),  ARLIM_P(phi_hi),
        Real* tmpr, ARLIM_P(tmpr_lo), ARLIM_P(tmpr_hi),
        const Real *h,
        const Real* tmpw, ARLIM_P(tmpw_lo), ARLIM_P(tmpw_hi),
        const Real* tmpp, ARLIM_P(tmpp_lo), ARLIM_P(tmpp_hi),
        const int* lo, const int* hi,
        const int* nc
        );

    void FORT_CGADVCP (
        Real* tmpp,       ARLIM_P(tmpp_lo), ARLIM_P(tmpp_hi),
        const Real* tmpr, ARLIM_P(tmpr_lo), ARLIM_P(tmpr_hi),
        const Real *b,
        const int* lo, const int* hi,
        const int* nc
        );

    void FORT_CGXDOTY (
        Real *pw,
        const Real* phi,  ARLIM_P(phi_lo),  ARLIM_P(phi_hi),
        const Real* tmpw, ARLIM_P(tmpw_lo), ARLIM_P(tmpw_hi),
        const int* lo, const int* hi,
        const int* nc
        );

  void FORT_CGSXAY (
        Real* s,  ARLIM_P(s_lo),  ARLIM_P(s_hi),
        const Real* x,  ARLIM_P(x_lo),  ARLIM_P(x_hi),
	const Real* a,
        const Real* y, ARLIM_P(y_lo), ARLIM_P(y_hi),
        const int* lo, const int* hi,
        const int* nc
        );

}
#endif

#endif /*_CG_F_H_*/
