#ifndef _DivVis_F_H_

#define _DivVis_F_H_ "$Id"

#include <REAL.H>


#if   defined(BL_LANG_FORT)
#   if  (BL_SPACEDIM==2)
#     define FORT_DVAPPLY          dvapply2d
#     define FORT_DVFLUX           dvflux2d
#     define FORT_GSRB             gsrbvisc2d
#     define FORT_APPLYBC          mcapplybc2dgen
#     define FORT_RESIDL           mcresid2dgen
#     define FORT_AVERAGECC        mcaveragecc2dgen
#     define FORT_AVERAGEEC        mcaverageec2dgen
#     define FORT_HARMONIC_AVERAGEEC mcharaverageec2dgen
#   elif(BL_SPACEDIM==3)
#     define FORT_DVAPPLY          dvapply3d
#     define FORT_DVFLUX           dvflux3d
#     define FORT_GSRB             gsrbvisc3d
#     define FORT_APPLYBC          mcapplybc3dgen
#     define FORT_RESIDL           mcresid3dgen
#     define FORT_AVERAGECC        mcaveragecc3dgen
#     define FORT_AVERAGEEC        mcaverageec3dgen
#     define FORT_HARMONIC_AVERAGEEC mcharaverageec3dgen
#   endif
#endif

#if   defined(BL_LANG_CC)

#if (BL_SPACEDIM == 2)

# ifdef BL_FORT_USE_UPPERCASE
#  define FORT_DVAPPLY          DVAPPLY2D
#  define FORT_DVFLUX           DVFLUX2D
#  define FORT_GSRB             GSRBVISC2D
#  define FORT_APPLYBC          MCAPPLYBC2DGEN
#  define FORT_RESIDL           MCRESID2DGEN
#  define FORT_AVERAGECC        MCAVERAGECC2DGEN
#  define FORT_AVERAGEEC        MCAVERAGEEC2DGEN
#  define FORT_HARMONIC_AVERAGEEC MCHARAVERAGEEC2DGEN
# endif

# ifdef BL_FORT_USE_UNDERSCORE
#  define FORT_DVAPPLY          dvapply2d_
#  define FORT_DVFLUX           dvflux2d_
#  define FORT_GSRB             gsrbvisc2d_
#  define FORT_APPLYBC          mcapplybc2dgen_
#  define FORT_RESIDL           mcresid2dgen_
#  define FORT_AVERAGECC        mcaveragecc2dgen_
#  define FORT_AVERAGEEC        mcaverageec2dgen_
#  define FORT_HARMONIC_AVERAGEEC mcharaverageec2dgen_
# endif

#else

# ifdef BL_FORT_USE_UPPERCASE
#  define FORT_DVAPPLY          DVAPPLY3D
#  define FORT_DVFLUX           DVFLUX3D
#  define FORT_GSRB             GSRBVISC3D
#  define FORT_APPLYBC          MCAPPLYBC3DGEN
#  define FORT_RESIDL           MCRESID3DGEN
#  define FORT_AVERAGECC        MCAVERAGECC3DGEN
#  define FORT_AVERAGEEC        MCAVERAGEEC3DGEN
#  define FORT_HARMONIC_AVERAGEEC MCHARAVERAGEEC3DGEN
# endif

# ifdef BL_FORT_USE_UNDERSCORE
#  define FORT_DVAPPLY          dvapply3d_
#  define FORT_DVFLUX           dvflux3d_
#  define FORT_GSRB             gsrbvisc3d_
#  define FORT_APPLYBC          mcapplybc3dgen_
#  define FORT_RESIDL           mcresid3dgen_
#  define FORT_AVERAGECC        mcaveragecc3dgen_
#  define FORT_AVERAGEEC        mcaverageec3dgen_
#  define FORT_HARMONIC_AVERAGEEC mcharaverageec3dgen_
# endif

#endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_AVERAGECC (
	Real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
	const Real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
	const int *tlo, const int *thi, const int *nc
	);

    void FORT_RESIDL (
	Real* res      , ARLIM_P(res_lo), ARLIM_P(res_hi),
	const Real* rhs, ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
	const Real* phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
	const int* lo, const int* hi, const int* nc
	);
    void FORT_AVERAGEEC (
	Real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
	const Real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
	const int *tlo, const int *thi, const int *nc,
	const int *axis
	);

    void FORT_HARMONIC_AVERAGEEC (
	Real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
	const Real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
	const int *tlo, const int *thi, const int *nc,
	const int *axis
	);
  void FORT_DVAPPLY(
		    const Real*u, ARLIM_P(ulo), ARLIM_P(uhi),
		    const Real*alpha, const Real*beta,
		    const Real*a, ARLIM_P(alo),   ARLIM_P(ahi),
		    const Real*muX, ARLIM_P(muXlo), ARLIM_P(muXhi),
		    const Real*muY, ARLIM_P(muYlo), ARLIM_P(muYhi),
#if (BL_SPACEDIM > 2)
		    const Real*muZ, ARLIM_P(muZlo), ARLIM_P(muZhi),
#endif
		    Real*out, ARLIM_P(outlo), ARLIM_P(outhi),
		    const int*maskn, ARLIM_P(masknlo), ARLIM_P(masknhi),
		    const int*maske, ARLIM_P(maskelo), ARLIM_P(maskehi),
		    const int*maskw, ARLIM_P(maskwlo), ARLIM_P(maskwhi),
		    const int*masks, ARLIM_P(maskslo), ARLIM_P(maskshi),
#if (BL_SPACEDIM > 2)
		    const int*maskt, ARLIM_P(masktlo), ARLIM_P(maskthi),
		    const int*maskb, ARLIM_P(maskblo), ARLIM_P(maskbhi),
#endif
		    const Real*trandern,
		    ARLIM_P(trandernlo),ARLIM_P(trandernhi),
		    const Real*trandere,
		    ARLIM_P(tranderelo),ARLIM_P(tranderehi),
		    const Real*tranderw,
		    ARLIM_P(tranderwlo),ARLIM_P(tranderwhi),
		    const Real*tranders,
		    ARLIM_P(tranderslo),ARLIM_P(trandershi),
#if (BL_SPACEDIM > 2)
		    const Real*trandert,
		    ARLIM_P(trandertlo),ARLIM_P(tranderthi),
		    const Real*tranderb,
		    ARLIM_P(tranderblo),ARLIM_P(tranderbhi),
#endif
		    const int* lo, const int* hi, const Real* h);

  void FORT_DVFLUX(
		   const Real*u, ARLIM_P(ulo), ARLIM_P(uhi),
		   const Real*alpha, const Real*beta,
		   const Real*a, ARLIM_P(alo),   ARLIM_P(ahi),
		   const Real*muX, ARLIM_P(muXlo), ARLIM_P(muXhi),
		   const Real*muY, ARLIM_P(muYlo), ARLIM_P(muYhi),
#if BL_SPACEDIM>2
		   const Real*muZ, ARLIM_P(muZlo), ARLIM_P(muZhi),
#endif
		   Real*xflux, ARLIM_P(xfluxlo), ARLIM_P(xfluxhi),
		   Real*yflux, ARLIM_P(yfluxlo), ARLIM_P(yfluxhi),
#if BL_SPACEDIM>2
		   Real*zflux, ARLIM_P(zfluxlo), ARLIM_P(zfluxhi),
#endif
		   const int*maskn, ARLIM_P(masknlo), ARLIM_P(masknhi),
		   const int*maske, ARLIM_P(maskelo), ARLIM_P(maskehi),
		   const int*maskw, ARLIM_P(maskwlo), ARLIM_P(maskwhi),
		   const int*masks, ARLIM_P(maskslo), ARLIM_P(maskshi),
#if BL_SPACEDIM>2
		   const int*maskt, ARLIM_P(masktlo), ARLIM_P(maskthi),
		   const int*maskb, ARLIM_P(maskblo), ARLIM_P(maskbhi),
#endif
		   const Real*trandern,
		   ARLIM_P(trandernlo),ARLIM_P(trandernhi),
		   const Real*trandere,
		   ARLIM_P(tranderelo),ARLIM_P(tranderehi),
		   const Real*tranderw,
		   ARLIM_P(tranderwlo),ARLIM_P(tranderwhi),
		   const Real*tranders,
		   ARLIM_P(tranderslo),ARLIM_P(trandershi),
#if BL_SPACEDIM>2
		   const Real*trandert,
		   ARLIM_P(trandertlo),ARLIM_P(tranderthi),
		   const Real*tranderb,
		   ARLIM_P(tranderblo),ARLIM_P(tranderbhi),
#endif
		   const int* lo, const int* hi, const Real* h);

  void FORT_GSRB(
		 Real* u, ARLIM_P(ulo), ARLIM_P(uhi),
		 const Real* rhs, ARLIM_P(rhslo), ARLIM_P(rhshi),
		 Real* alpha, Real* beta,
		 const Real* a, ARLIM_P(alo), ARLIM_P(ahi),
		 const Real* muX, ARLIM_P(muXlo), ARLIM_P(muXhi),
		 const Real* muY, ARLIM_P(muYlo), ARLIM_P(muYhi),
#if BL_SPACEDIM>2
		 const Real* muZ, ARLIM_P(muZlo), ARLIM_P(muZhi),
#endif
		 const int* maskn, ARLIM_P(masknlo), ARLIM_P(masknhi),
		 const Real* fn, ARLIM_P(fnlo), ARLIM_P(fnhi),
		 const int* maske, ARLIM_P(maskelo), ARLIM_P(maskehi),
		 const Real* fe, ARLIM_P(felo), ARLIM_P(fehi),
		 const int* maskw, ARLIM_P(maskwlo), ARLIM_P(maskwhi),
		 const Real* fw, ARLIM_P(fwlo), ARLIM_P(fwhi),
		 const int* masks, ARLIM_P(maskslo), ARLIM_P(maskshi),
		 const Real* fs, ARLIM_P(fslo), ARLIM_P(fshi),
#if BL_SPACEDIM>2
		 const int* maskt, ARLIM_P(masktlo), ARLIM_P(maskthi),
		 const Real* ft, ARLIM_P(ftlo), ARLIM_P(fthi),
		 const int* maskb, ARLIM_P(maskblo), ARLIM_P(maskbhi),
		 const Real* fb, ARLIM_P(fblo), ARLIM_P(fbhi),
#endif
		 const Real* trandern, 
		 ARLIM_P(trandernlo), ARLIM_P(trandernhi),
		 const Real* trandere, 
		 ARLIM_P(tranderelo), ARLIM_P(tranderehi),
		 const Real* tranderw, 
		 ARLIM_P(tranderwlo), ARLIM_P(tranderwhi),
		 const Real* tranders, 
		 ARLIM_P(tranderslo), ARLIM_P(trandershi),
#if BL_SPACEDIM>2
		 const Real* trandert, 
		 ARLIM_P(trandertlo), ARLIM_P(tranderthi),
		 const Real* tranderb, 
		 ARLIM_P(tranderblo), ARLIM_P(tranderbhi),
#endif
		 const int* lo, const int* hi,
		 const Real* h, int &nc, int &phaseflag);
		 
  void FORT_APPLYBC(
		    const int *flagden,  // 1 if want values in den
		    const int *flagbc,   // 1 for inhomogeneous
		    const int *maxorder, // maxorder of interpolant
		    Real *phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
		    const int *cdir,
		    const int *bct,
		    const Real *bcl,
		    const Real *bcval, ARLIM_P(bcval_lo), ARLIM_P(bcval_hi),
#if BL_SPACEDIM==2
		    const int *mask,   ARLIM_P(mask_lo),  ARLIM_P(mask_hi), 
		    const int *maskp,   ARLIM_P(maskp_lo),  ARLIM_P(maskp_hi), 
		    const int *maskm,   ARLIM_P(maskm_lo),  ARLIM_P(maskm_hi), 
#elif BL_SPACEDIM==3
		    const int *maskn,   ARLIM_P(maskn_lo),  ARLIM_P(maskn_hi), 
		    const int *maske,   ARLIM_P(maske_lo),  ARLIM_P(maske_hi), 
		    const int *maskw,   ARLIM_P(maskw_lo),  ARLIM_P(maskw_hi), 
		    const int *masks,   ARLIM_P(masks_lo),  ARLIM_P(masks_hi), 
		    const int *maskt,   ARLIM_P(maskt_lo),  ARLIM_P(maskt_hi), 
		    const int *maskb,   ARLIM_P(maskb_lo),  ARLIM_P(maskb_hi), 
#endif
		    Real *den,         ARLIM_P(den_lo),   ARLIM_P(den_hi),
		    const Real *exttd, ARLIM_P(exttd_lo), ARLIM_P(exttd_hi),
		    Real *tander,      ARLIM_P(tander_lo),ARLIM_P(tander_hi),
		    const int *lo, const int *hi, 
		    const int *nc,
		    const Real *h
		    );


};

#endif

#endif
