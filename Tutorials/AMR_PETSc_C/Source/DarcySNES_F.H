#ifndef _DarcySNES_F_H
#define _DarcySNES_F_H

#ifdef BL_LANG_FORT
#    define FORT_RS_PDOTRES  rspdotres
#    define FORT_DARCYFLUX   rsdarcyflux
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_RS_PDOTRES  RSPDOTRES
#    define FORT_DARCYFLUX   RSDARCYFLUX
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_RS_PDOTRES  rspdotres
#    define FORT_DARCYFLUX   rsdarcyflux
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_RS_PDOTRES  rspdotres_
#    define FORT_DARCYFLUX   rsdarcyflux_
#  endif

#include <SPACE.H>
#include <ArrayLim.H>

extern "C"
{
  void FORT_RS_PDOTRES (Real* res,        ARLIM_P(reslo), ARLIM_P(reshi),
			const Real* rso,  ARLIM_P(rsnlo), ARLIM_P(rsnhi),
			const Real* rsn,  ARLIM_P(rsnlo), ARLIM_P(rsnhi),
			const Real* dt, const Real* phi,
                        const int*  lo, const int*  hi, const int* nc);

  void FORT_DARCYFLUX(const Real* p,  ARLIM_P(plo),  ARLIM_P(phi),
		      const Real* kr, ARLIM_P(krlo), ARLIM_P(krhi),
		      const int*  i,  ARLIM_P(ilo),  ARLIM_P(ihi),
		      D_DECL(Real* fx, Real* fy, Real* fz),
		      D_DECL(ARLIM_P(fxlo), ARLIM_P(fylo), ARLIM_P(fzlo)),
		      D_DECL(ARLIM_P(fxhi), ARLIM_P(fyhi), ARLIM_P(fzhi)),
		      const int* lo, const int* hi, const Real* kappax, const Real* kappaz,
		      const Real* rho, const Real* mu, const Real* gravity, const Real* dx);
}

#endif
#endif /*_RSARDSOLVER_F_H*/
