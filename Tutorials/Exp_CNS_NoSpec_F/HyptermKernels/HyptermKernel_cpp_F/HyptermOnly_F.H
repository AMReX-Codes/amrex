/* -------------------------------------------------------------------- */
/*  HyptermOnly_F.H                                                     */
/* -------------------------------------------------------------------- */
#ifndef _HYPTERMONLY_F_H_
#define _HYPTERMONLY_F_H_

#ifdef BL_LANG_FORT
#    define FORT_HYPTERM     hypterm
#    define FORT_HYPTERM_UNOPT     hypterm_unopt
#    define FORT_INITDATA     initdata
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_HYPTERM     HYPTERM
#    define FORT_HYPTERM_UNOPT     HYPTERM_UNOPT
#    define FORT_INITDATA     INITDATA
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_HYPTERM     hypterm
#    define FORT_HYPTERM_UNOPT     hypterm_unopt
#    define FORT_INITDATA     initdata
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_HYPTERM     hypterm_
#    define FORT_HYPTERM_UNOPT     hypterm_unopt_
#    define FORT_INITDATA     initdata_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_HYPTERM(const Real *cons, ARLIM_P(cons_lo), ARLIM_P(cons_hi),
		      ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                      const Real *dx,
                      const Real *qx, const Real *flux,
                      const int *nComp);

    void FORT_HYPTERM_UNOPT(const Real *cons, ARLIM_P(cons_lo), ARLIM_P(cons_hi),
		      ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                      const Real *dx,
                      const Real *qx, const Real *flux,
                      const int *nComp);

    void FORT_INITDATA(const Real *cons, ARLIM_P(cons_lo), ARLIM_P(cons_hi),
                       const Real *q, const Real *dx, const int *nComp);
}
#endif

#endif
