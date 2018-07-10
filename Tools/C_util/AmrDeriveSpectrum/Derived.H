
#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>

#ifdef BL_LANG_FORT
#  define   FORT_DERLGDN       fdlgdn3d
#  define   FORT_DERDEN        fdden3d
#  define   FORT_DERTRAC       fdtrac3d
#  define   FORT_DERTMASK      fdtmask3d
#  define   FORT_DERSTATE      fdstate3d
#  define   FORT_DERXVEL       fdxvel3d
#  define   FORT_DERYVEL       fdyvel3d
#  define   FORT_DERZVEL       fdzvel3d
#  define   FORT_DERTENG       fdeng3d
#  define   FORT_DERIENG       fdieng3d
#  define   FORT_DERKENG       fdkeng3d
#  define   FORT_DERVORT       fdvort3d
#  define   FORT_DERVORTCOMPS       fdvortcomps3d
#  define   FORT_DERCON        fdcon3d
#  define   FORT_INTEGRATE     fdintegrate
#  define   FORT_JPDF          fdjpdf
#  define   FORT_LIMITS        fdlimits
#  define   FORT_SETGAMMA      fdsetgam
#  define   FORT_SETALPHA      fdsetalpha
#  define   FORT_SETQPOS       fdsetqpos
#  define   FORT_SUMMASS       fdsummass
#  define   FORT_SUMCONC       fdsumconc
#  define   FORT_SUMVORT       fdsumvort
#  define   FORT_DERDILI       fddili3d
#  define   FORT_BINTRAC       bintrac
#  define   FORT_DERMACH       fdmach3d
#  define   FORT_RANDMOVE      randmove
#  define   FORT_CHEMMOVE      chemmove
#  define   FORT_MAXLAMBDA     maxlambda
#  define   FORT_MAXDTLAMBDA   maxdtlambda
#  define   FORT_MAXDTDIFF     maxdtdiff
#  define   FORT_SELECTPOS     selectpos
#  define   FORT_PROBFAB       probfab
#  define   FORT_ACCUMPROB     accumprob
#  define   FORT_TEMPQMODEL    tempqmodel
#  define   FORT_FUELFRAC      fuelfrac
#  define   FORT_OXFRAC        oxfrac
#  define   FORT_PRODFRAC      prodfrac
#else

#  if defined(BL_FORT_USE_UPPERCASE)
#    define   FORT_DERLGDN       FDLGDN3D
#    define   FORT_DERDEN        FDDEN3D
#    define   FORT_DERTRAC       FDTRAC3D
#    define   FORT_DERTMASK      FDTMASK3D
#    define   FORT_DERSTATE      FDSTATE3D
#    define   FORT_DERXVEL       FDXVEL3D
#    define   FORT_DERYVEL       FDYVEL3D
#    define   FORT_DERZVEL       FDZVEL3D
#    define   FORT_DERTENG       FDENG3D
#    define   FORT_DERIENG       FDIENG3D
#    define   FORT_DERKENG       FDKENG3D
#    define   FORT_DERVORT       FDVORT3D
#    define   FORT_DERVORTCOMPS       FDVORTCOMPS3D
#    define   FORT_DERCON        FDCON3D
#    define   FORT_INTEGRATE     FDINTEGRATE
#    define   FORT_JPDF          FDJPDF
#    define   FORT_LIMITS        FDLIMITS
#    define   FORT_SETGAMMA      FDSETGAM
#    define   FORT_SETALPHA      FDSETALPHA
#    define   FORT_SETQPOS       FDSETQPOS
#    define   FORT_SUMMASS       FDSUMMASS
#    define   FORT_SUMCONC       FDSUMCONC
#    define   FORT_SUMVORT       FDSUMVORT
#    define   FORT_DERDILI       FDDILI3D
#    define   FORT_BINTRAC       BINTRAC
#    define   FORT_DERMACH       FDMACH3D
#    define   FORT_RANDMOVE      RANDMOVE
#    define   FORT_CHEMMOVE      CHEMMOVE
#    define   FORT_MAXLAMBDA     MAXLAMBDA
#    define   FORT_MAXDTLAMBDA   MAXDTLAMBDA
#    define   FORT_MAXDTDIFF     MAXDTDIFF
#    define   FORT_SELECTPOS     SELECTPOS
#    define   FORT_BUILDPROB     BUILDPROB
#    define   FORT_PROBFAB       PROBFAB
#    define   FORT_ACCUMPROB     ACCUMPROB
#    define   FORT_TEMPQMODEL    TEMPQMODEL
#    define   FORT_FUELFRAC      FUELFRAC
#    define   FORT_OXFRAC        OXFRAC
#    define   FORT_PRODFRAC      PRODFRAC
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define   FORT_DERLGDN       fdlgdn3d
#    define   FORT_DERDEN        fdden3d
#    define   FORT_DERTRAC       fdtrac3d
#    define   FORT_DERTMASK      fdtmask3d
#    define   FORT_DERSTATE      fdstate3d
#    define   FORT_DERXVEL       fdxvel3d
#    define   FORT_DERYVEL       fdyvel3d
#    define   FORT_DERZVEL       fdzvel3d
#    define   FORT_DERTENG       fdeng3d
#    define   FORT_DERIENG       fdieng3d
#    define   FORT_DERKENG       fdkeng3d
#    define   FORT_DERVORT       fdvort3d
#    define   FORT_DERVORTCOMPS       fdvortcomps3d
#    define   FORT_DERCON        fdcon3d
#    define   FORT_INTEGRATE     fdintegrate
#    define   FORT_JPDF          fdjpdf
#    define   FORT_LIMITS        fdlimits
#    define   FORT_SETGAMMA      fdsetgam
#    define   FORT_SETALPHA      fdsetalpha
#    define   FORT_SETQPOS       fdsetqpos
#    define   FORT_SUMMASS       fdsummass
#    define   FORT_SUMCONC       fdsumconc
#    define   FORT_SUMVORT       fdsumvort
#    define   FORT_DERDILI       fddili3d
#    define   FORT_BINTRAC       bintrac
#    define   FORT_DERMACH       fdmach3d
#    define   FORT_DERMACH       fdmach3d
#    define   FORT_RANDMOVE      randmove
#    define   FORT_CHEMMOVE      chemmove
#    define   FORT_MAXLAMBDA     maxlambda
#    define   FORT_MAXDTLAMBDA   maxdtlambda
#    define   FORT_MAXDTDIFF     maxdtdiff
#    define   FORT_SELECTPOS     selectpos
#    define   FORT_PROBFAB       probfab
#    define   FORT_ACCUMPROB     accumprob
#    define   FORT_TEMPQMODEL    tempqmodel
#    define   FORT_FUELFRAC      fuelfrac
#    define   FORT_OXFRAC        oxfrac
#    define   FORT_PRODFRAC      prodfrac
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define   FORT_DERLGDN       fdlgdn3d_
#    define   FORT_DERDEN        fdden3d_
#    define   FORT_DERTRAC       fdtrac3d_
#    define   FORT_DERTMASK      fdtmask3d_
#    define   FORT_DERSTATE      fdstate3d_
#    define   FORT_DERXVEL       fdxvel3d_
#    define   FORT_DERYVEL       fdyvel3d_
#    define   FORT_DERZVEL       fdzvel3d_
#    define   FORT_DERTENG       fdeng3d_
#    define   FORT_DERIENG       fdieng3d_
#    define   FORT_DERKENG       fdkeng3d_
#    define   FORT_DERVORT       fdvort3d_
#    define   FORT_DERVORTCOMPS       fdvortcomps3d_
#    define   FORT_DERCON        fdcon3d_
#    define   FORT_INTEGRATE     fdintegrate_
#    define   FORT_JPDF          fdjpdf_
#    define   FORT_LIMITS        fdlimits_
#    define   FORT_SETGAMMA      fdsetgam_
#    define   FORT_SETALPHA      fdsetalpha_
#    define   FORT_SETQPOS       fdsetqpos_
#    define   FORT_SUMMASS       fdsummass_
#    define   FORT_SUMCONC       fdsumconc_
#    define   FORT_SUMVORT       fdsumvort_
#    define   FORT_DERDILI       fddili3d_
#    define   FORT_BINTRAC       bintrac_
#    define   FORT_DERMACH       fdmach3d_
#    define   FORT_DERMACH       fdmach3d_
#    define   FORT_RANDMOVE      randmove_
#    define   FORT_CHEMMOVE      chemmove_
#    define   FORT_MAXLAMBDA     maxlambda_
#    define   FORT_MAXDTLAMBDA   maxdtlambda_
#    define   FORT_MAXDTDIFF     maxdtdiff_
#    define   FORT_SELECTPOS     selectpos_
#    define   FORT_PROBFAB       probfab_
#    define   FORT_ACCUMPROB     accumprob_
#    define   FORT_TEMPQMODEL    tempqmodel_
#    define   FORT_FUELFRAC      fuelfrac_
#    define   FORT_OXFRAC        oxfrac_
#    define   FORT_PRODFRAC      prodfrac_
#  endif

//
// Use ARLIM in the following.
//

extern "C"
{
  void FORT_SETGAMMA( const Real *gam );

  void FORT_SETALPHA( const Real *alp );

  void FORT_SETQPOS( const int *qpos );

  void FORT_DERLGDN(const Real *u, 
		    const int &uloi1, const int &uloi2, const int &uloi3,
		    const int &uhii1, const int &uhii2, const int &uhii3,
		    const int &nvar,
		    const Real *xlo, const Real *xhi,
		    Real *dat,
		    const int &dloi1, const int &dloi2, const int &dloi3,
		    const int &dhii1, const int &dhii2, const int &dhii3,
		    const int *ovlo, const int *ovhi);

  void FORT_DERXVEL(const Real *u, 
		    const int &uloi1, const int &uloi2, const int &uloi3,
		    const int &uhii1, const int &uhii2, const int &uhii3,
		    const int &nvar,
		    const Real *xlo, const Real *xhi,
		    Real *dat,
		    const int &dloi1, const int &dloi2, const int &dloi3,
		    const int &dhii1, const int &dhii2, const int &dhii3,
		    const int *ovlo, const int *ovhi);

  void FORT_DERYVEL(const Real *u, 
		    const int &uloi1, const int &uloi2, const int &uloi3,
		    const int &uhii1, const int &uhii2, const int &uhii3,
		    const int &nvar,
		    const Real *xlo, const Real *xhi,
		    Real *dat,
		    const int &dloi1, const int &dloi2, const int &dloi3,
		    const int &dhii1, const int &dhii2, const int &dhii3,
		    const int *ovlo, const int *ovhi);

  void FORT_DERZVEL(const Real *u, 
		    const int &uloi1, const int &uloi2, const int &uloi3,
		    const int &uhii1, const int &uhii2, const int &uhii3,
		    const int &nvar,
		    const Real *xlo, const Real *xhi,
		    Real *dat,
		    const int &dloi1, const int &dloi2, const int &dloi3,
		    const int &dhii1, const int &dhii2, const int &dhii3,
		    const int *ovlo, const int *ovhi);

  void FORT_DERIENG(const Real *u,
		    const int &uloi1, const int &uloi2, const int &uloi3,
		    const int &uhii1, const int &uhii2, const int &uhii3,
		    const int &nvar,
		    const Real *xlo, const Real *xhi,
		    Real *dat,
		    const int &dloi1, const int &dloi2, const int &dloi3,
		    const int &dhii1, const int &dhii2, const int &dhii3,
		    const int *ovlo, const int *ovhi);

  void FORT_DERCON(const Real *u, 
		   const int &uloi1, const int &uloi2, const int &uloi3,
		   const int &uhii1, const int &uhii2, const int &uhii3,
		   const int &nvar,
		   const Real *xlo, const Real *xhi,
		   Real *dat,
		   const int &dloi1, const int &dloi2,const int &dloi3,
		   const int &dhii1, const int &dhii2,const int &dhii3,
		   const int *ovlo, const int *ovhi);

  void FORT_DERVORT(const Real *u, 
                    const int &uloi1, const int &uloi2, const int &uloi3,
                    const int &uhii1, const int &uhii2, const int &uhii3,
                    const int &nvar,
                    const Real *xlo, const Real *xhi,
                    Real *dat,
                    const int &dloi1, const int &dloi2,const int &dloi3,
                    const int &dhii1, const int &dhii2,const int &dhii3,
                    const int *ovlo, const int *ovhi);

  void FORT_DERVORTCOMPS(const Real *u, 
                    const int &uloi1, const int &uloi2, const int &uloi3,
                    const int &uhii1, const int &uhii2, const int &uhii3,
                    const int &nvar,
                    const Real *xlo, const Real *xhi,
                    Real *dat,
                    const int &dloi1, const int &dloi2,const int &dloi3,
                    const int &dhii1, const int &dhii2,const int &dhii3,
                    const int *ovlo, const int *ovhi);

  void FORT_TEMPQMODEL(const Real *u, 
                    const int &uloi1, const int &uloi2, const int &uloi3,
                    const int &uhii1, const int &uhii2, const int &uhii3,
                    const int &nvar,
                    const Real *xlo, const Real *xhi,
                    Real *dat,
                    const int &dloi1, const int &dloi2,const int &dloi3,
                    const int &dhii1, const int &dhii2,const int &dhii3,
                    const int *ovlo, const int *ovhi);

  void FORT_FUELFRAC(const Real *u, 
                    const int &uloi1, const int &uloi2, const int &uloi3,
                    const int &uhii1, const int &uhii2, const int &uhii3,
                    const int &nvar,
                    const Real *xlo, const Real *xhi,
                    Real *dat,
                    const int &dloi1, const int &dloi2,const int &dloi3,
                    const int &dhii1, const int &dhii2,const int &dhii3,
                    const int *ovlo, const int *ovhi);

  void FORT_OXFRAC(const Real *u, 
                    const int &uloi1, const int &uloi2, const int &uloi3,
                    const int &uhii1, const int &uhii2, const int &uhii3,
                    const int &nvar,
                    const Real *xlo, const Real *xhi,
                    Real *dat,
                    const int &dloi1, const int &dloi2,const int &dloi3,
                    const int &dhii1, const int &dhii2,const int &dhii3,
                    const int *ovlo, const int *ovhi);

  void FORT_PRODFRAC(const Real *u, 
                    const int &uloi1, const int &uloi2, const int &uloi3,
                    const int &uhii1, const int &uhii2, const int &uhii3,
                    const int &nvar,
                    const Real *xlo, const Real *xhi,
                    Real *dat,
                    const int &dloi1, const int &dloi2,const int &dloi3,
                    const int &dhii1, const int &dhii2,const int &dhii3,
                    const int *ovlo, const int *ovhi);

  void FORT_DERKENG(const Real *u, 
                    const int &uloi1, const int &uloi2, const int &uloi3,
                    const int &uhii1, const int &uhii2, const int &uhii3,
                    const int &nvar,
                    const Real *xlo, const Real *xhi,
                    Real *dat,
                    const int &dloi1, const int &dloi2,const int &dloi3,
                    const int &dhii1, const int &dhii2,const int &dhii3,
                    const int *ovlo, const int *ovhi);

  void FORT_DERTENG(const Real *u, 
                    const int &uloi1, const int &uloi2, const int &uloi3,
                    const int &uhii1, const int &uhii2, const int &uhii3,
                    const int &nvar,
                    const Real *xlo, const Real *xhi,
                    Real *dat,
                    const int &dloi1, const int &dloi2,const int &dloi3,
                    const int &dhii1, const int &dhii2,const int &dhii3,
                    const int *ovlo, const int *ovhi);

  void FORT_DERDILI(const Real *u, 
                    const int &uloi1, const int &uloi2, const int &uloi3,
                    const int &uhii1, const int &uhii2, const int &uhii3,
                    const int &nvar,
                    const Real *xlo, const Real *xhi,
                    Real *dat,
                    const int &dloi1, const int &dloi2,const int &dloi3,
                    const int &dhii1, const int &dhii2,const int &dhii3,
                    const int *ovlo, const int *ovhi);

   void FORT_INTEGRATE (const Real* field, 
			 ARLIM_P(field_lo), ARLIM_P(field_hi),
			 ARLIM_P(valid_lo), ARLIM_P(valid_hi),
			 const Real* delta, 
			 Real* integral);

    void FORT_JPDF (const Real* field, 
		    const Real* visible,
		    ARLIM_P(field_lo), ARLIM_P(field_hi),
		    ARLIM_P(valid_lo), ARLIM_P(valid_hi),
		    const int* hBins, const Real* hMax, const Real* hMin,
		    const int* vBins, const Real* vMax, const Real* vMin,
		    int* excessCount, int* jointCount);

    void FORT_LIMITS (const Real* field, 
		      const Real* visible,
		      ARLIM_P(field_lo), ARLIM_P(field_hi),
		      ARLIM_P(valid_lo), ARLIM_P(valid_hi),
		      Real* maximum,
		      Real* minimum);

#if (BL_SPACEDIM == 2)
    void FORT_SUMMASS (const Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const Real* delta, Real* mass,
                       const Real* radius,
                       const int* irlo, const int* irhi,
                       const int* rz_flag, Real* tmp);
#else
    void FORT_SUMMASS (const Real* dat, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
                       ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                       const Real* delta, Real* mass, Real* tmp);
#endif

  void FORT_SUMCONC (const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
                     ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                     const Real* delta, Real* mass, Real* tmp);

  void FORT_SUMVORT (const Real* dat, ARLIM_P(dlo), ARLIM_P(dhi),
                     ARLIM_P(grid_lo), ARLIM_P(grid_hi),
                     const Real* delta, Real* mass, Real* tmp);

  void FORT_BINTRAC (const Real* trac, ARLIM_P(trac_lo), ARLIM_P(trac_hi),
                     const Real* rho,
                     const int* lo, const int* hi,
                     const Real* delta, Real* sum);

  void FORT_DERDEN (const Real *u, 
		    const int &uloi1, const int &uloi2, const int &uloi3,
		    const int &uhii1, const int &uhii2, const int &uhii3,
		    const int &nvar,
		    const Real *xlo, const Real *xhi,
		    Real *dat,
		    const int &dloi1, const int &dloi2, const int &dloi3,
		    const int &dhii1, const int &dhii2, const int &dhii3,
		    const int *ovlo, const int *ovhi);

  void FORT_DERTRAC(const Real *u, 
		    const int &uloi1, const int &uloi2, const int &uloi3,
		    const int &uhii1, const int &uhii2, const int &uhii3,
		    const int &nvar,
		    const Real *xlo, const Real *xhi,
		    Real *dat,
		    const int &dloi1, const int &dloi2, const int &dloi3,
		    const int &dhii1, const int &dhii2, const int &dhii3,
		    const int *ovlo, const int *ovhi);

  void FORT_DERTMASK(const Real *u, 
		    const int &uloi1, const int &uloi2, const int &uloi3,
		    const int &uhii1, const int &uhii2, const int &uhii3,
		    const int &nvar,
		    const Real *xlo, const Real *xhi,
		    Real *dat,
		    const int &dloi1, const int &dloi2, const int &dloi3,
		    const int &dhii1, const int &dhii2, const int &dhii3,
		    const int *ovlo, const int *ovhi);

  void FORT_DERSTATE(const Real *u, 
		    const int &uloi1, const int &uloi2, const int &uloi3,
		    const int &uhii1, const int &uhii2, const int &uhii3,
		    const int &nvar,
		    const Real *xlo, const Real *xhi,
		    Real *dat,
		    const int &dloi1, const int &dloi2, const int &dloi3,
		    const int &dhii1, const int &dhii2, const int &dhii3,
		    const int *ovlo, const int *ovhi);

  void FORT_DERMACH (const Real *u, 
		    const int &uloi1, const int &uloi2, const int &uloi3,
		    const int &uhii1, const int &uhii2, const int &uhii3,
		    const int &nvar,
		    const Real *xlo, const Real *xhi,
		    Real *dat,
		    const int &dloi1, const int &dloi2, const int &dloi3,
		    const int &dhii1, const int &dhii2, const int &dhii3,
		    const int *ovlo, const int *ovhi);

  void FORT_RANDMOVE (const Real* A,ARLIM_P(A_lo),ARLIM_P(A_hi),
                      const Real* B,
                      const Real* C,
                      const Real* D,
                      const Real* U,
                      const int*  nspec,
                      const Real* x, const Real* y, const int* spec,
                      const Real* dx, const Real* dy, const Real* dt,
                      const Real* dxinv, const Real* dyinv);

  void FORT_CHEMMOVE (const Real* RF,ARLIM_P(rf_lo),ARLIM_P(rf_hi),
                      const Real* RR,
                      const Real* CONC,
                      const int*  nspec,
                      const int*  nreac,
                      const int*  edges, const int* edgeslen,
                      const int*  pedges,
                      const Real* x, const Real* y, const int* spec,
                      const Real* dx, const Real* dy, const Real* dt,
                      const int* rxn,const Real* dxinv, const Real* dyinv);

  void FORT_MAXLAMBDA (const Real* RF,ARLIM_P(rf_lo),ARLIM_P(rf_hi),
                       const Real* RR,
                       const Real* CONC,
                       const int*  nspec,
                       const int*  nreac,
                       const int*  edges, const int* edgeslen,
                       const int*  pedges);

  void FORT_MAXDTLAMBDA (const Real* RF,ARLIM_P(rf_lo),ARLIM_P(rf_hi),
                         const Real* RR,
                         const Real* CONC,
                         const int*  nspec,
                         const int*  nreac,
                         const int*  edges, const int* edgeslen,
                         const int*  pedges,
                         const Real* x, const Real* y, const int* spec,
                         const Real* dx, const Real* dy, const Real* dt,
                         const Real* dxinv, const Real* dyinv);

  void FORT_MAXDTDIFF (const Real* C,ARLIM_P(C_lo),ARLIM_P(C_hi),
                       const int*  nspec,
                       const Real* x, const Real* y, const int* spec,
                       const Real* dx, const Real* dy, const Real* dt,
                       const Real* dxinv, const Real* dyinv);

  void FORT_PROBFAB (const Real* rf, ARLIM_P(rf_lo), ARLIM_P(rf_hi),
                     const Real* rr,
                     const Real* prob, ARLIM_P(prob_lo), ARLIM_P(prob_hi),
                     const int*  nspec,
                     const int*  nreac,
                     const int*  ispec,
                     const int*  edges, const int* edgeslen,
                     const int*  pedges,
                     const int*  isrz,
                     D_DECL(const Real* dx, const Real* dy, const Real* dz));

  void FORT_ACCUMPROB (const Real* prob, ARLIM_P(prob_lo), ARLIM_P(prob_hi),
                       const Real* totreact, Real* cumprob);

  void FORT_SELECTPOS (const Real* prob,ARLIM_P(prob_lo),ARLIM_P(prob_hi),
                       const Real* rn,
                       D_DECL(const Real* dx, const Real* dy, const Real* dz),
                       D_DECL(const Real* x, const Real* y, const Real* z));
  //
  // Typedef needed for defining a pointer to the FORTRAN functions
  // used in deriving other thermodynamic quantities from DATABOX
  // state variables.
  //
  typedef void (*DERIVE_FUNC) (const Real *u, 
                               const int &uloi1, const int &uloi2, const int &uloi3,
                               const int &uhii1, const int &uhii2, const int &uhii3,
                               const int &nvar,
                               const Real *xlo, const Real *xhi,
                               Real *dat,
                               const int &dloi1, const int &dloi2, const int &dloi3,
                               const int &dhii1, const int &dhii2, const int &dhii3,
                               const int *ovlo, const int *ovhi);
}


//
// Number of variables in the state including derived quantities.
//

#ifdef BL_DERIVE_IAMR
#define NSTATECOMPS 5
#else
#ifdef BL_DERIVE_HCHEM
#define NSTATECOMPS 14
#else
#define NSTATECOMPS 6
#endif
#endif

struct DeriveRecord
{
      std::string     name;  // Name of quantity to derive from state data.
      DERIVE_FUNC func;  // Pointer to function that computes this quantity.
      int         derMap[NSTATECOMPS];  // Map of required variables for the
				        // derived functions
					// 1 == needed,  0 == not needed
};

//
// The state quantities.
//

#if BL_SPACEDIM==3

#ifdef BL_DERIVE_IAMR

// nstates = spacedim + 2
// x_velocity
// y_velocity
// z_velocity
// density
// tracer

extern const struct DeriveRecord DerivedQuants[] =
{
  //  { "logden",        FORT_DERLGDN, { 0, 0, 0, 1, 0 } },
  //  { "x_velocity",    FORT_DERXVEL, { 1, 0, 0, 0, 0 } },
  //  { "y_velocity",    FORT_DERYVEL, { 0, 1, 0, 0, 0 } },
  //  { "z_velocity",    FORT_DERZVEL, { 0, 0, 1, 0, 0 } },
    { "density",       FORT_DERDEN,  { 0, 0, 0, 1, 0 } },
    { "tracer",        FORT_DERTRAC, { 0, 0, 0, 0, 1 } }
};

/*
// these are for iamr:  from the plotfile header:
// NavierStokes-V1.1
// 5
// x_velocity
// y_velocity
// z_velocity
// density
// tracer

extern const struct DeriveRecord DerivedQuants[] =
{
    { "logden",        FORT_DERLGDN, { 0, 0, 0, 1, 0 } },
    { "x_velocity",    FORT_DERXVEL, { 1, 0, 0, 0, 0 } },
    { "y_velocity",    FORT_DERYVEL, { 0, 1, 0, 0, 0 } },
    { "z_velocity",    FORT_DERZVEL, { 0, 0, 1, 0, 0 } },
    //{ "teng",          FORT_DERTENG, { 1, 0, 0, 0, 1 } },
    //{ "eint",          FORT_DERIENG, { 1, 1, 1, 1, 1 } },
    //{ "kint",          FORT_DERKENG, { 1, 1, 1, 1, 0 } },
    //{ "vorticity",     FORT_DERVORT, { 1, 1, 1, 1, 0 } },
    //{ "concentration", FORT_DERCON,  { 1, 0, 0, 0, 0 } },
    //{ "dilitation",    FORT_DERDILI, { 1, 1, 1, 1, 0 } },
    { "density",       FORT_DERDEN,  { 0, 0, 0, 1, 0 } },
    { "tracer",        FORT_DERTRAC, { 0, 0, 0, 0, 1 } }
};
*/

#else



#ifdef BL_DERIVE_HCHEM

// density
// xmom
// ymom
// zmom
// eden
// tracer
// rho.Y(O2)
// rho.Y(H2O)
// rho.Y(CO2)
// rho.Y(C2H2)
// rho.Y(CO)
// rho.Y(N2)
// pressure

extern const struct DeriveRecord DerivedQuants[] =
{
    { "logden",        FORT_DERLGDN, { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { "x_velocity",    FORT_DERXVEL, { 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { "y_velocity",    FORT_DERYVEL, { 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { "z_velocity",    FORT_DERZVEL, { 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    //{ "teng",          FORT_DERTENG, { 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 } },
    //{ "eint",          FORT_DERIENG, { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 } },
    //{ "kint",          FORT_DERKENG, { 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    //{ "vorticity",     FORT_DERVORT, { 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    //{ "concentration", FORT_DERCON,  { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 } },
    //{ "dilitation",    FORT_DERDILI, { 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { "density",       FORT_DERDEN,  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    //{ "mach",          FORT_DERMACH, { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { "state",          FORT_DERSTATE, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { "tmask",         FORT_DERTMASK, { 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0 } },
    { "tempqmodel",    FORT_TEMPQMODEL, { 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0 } },
    { "fuelfrac",      FORT_FUELFRAC, { 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 } },
    { "oxfrac",        FORT_OXFRAC,   { 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 } },
    { "prodfrac",      FORT_PRODFRAC, { 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 } },
    { "tracer",        FORT_DERTRAC,  { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 } }
};

#else

/*
extern const struct DeriveRecord DerivedQuants[] =
{
    { "logden",        FORT_DERLGDN, { 1, 0, 0, 0, 0, 0 } },
    { "x_velocity",    FORT_DERXVEL, { 1, 1, 0, 0, 0, 0 } },
    { "y_velocity",    FORT_DERYVEL, { 1, 0, 1, 0, 0, 0 } },
    { "z_velocity",    FORT_DERZVEL, { 1, 0, 0, 1, 0, 0 } },
    { "teng",          FORT_DERTENG, { 1, 0, 0, 0, 1, 0 } },
    { "eint",          FORT_DERIENG, { 1, 1, 1, 1, 1, 0 } },
    { "kint",          FORT_DERKENG, { 1, 1, 1, 1, 0, 0 } },
    { "vorticity",     FORT_DERVORT, { 1, 1, 1, 1, 0, 0 } },
    { "concentration", FORT_DERCON,  { 1, 0, 0, 0, 0, 1 } },
    { "dilitation",    FORT_DERDILI, { 1, 1, 1, 1, 0, 0 } },
    { "density",       FORT_DERDEN,  { 1, 0, 0, 0, 0, 0 } },
    { "tracer",        FORT_DERTRAC, { 0, 0, 0, 0, 0, 1 } }
};
*/

extern const struct DeriveRecord DerivedQuants[] =
{

    { "logden",        FORT_DERLGDN, { 1, 0, 0, 0, 0, 0 } },
    { "x_velocity",    FORT_DERXVEL, { 1, 1, 0, 0, 0, 0 } },
    { "y_velocity",    FORT_DERYVEL, { 1, 0, 1, 0, 0, 0 } },
    { "z_velocity",    FORT_DERZVEL, { 1, 0, 0, 1, 0, 0 } },
    { "teng",          FORT_DERTENG, { 1, 0, 0, 0, 1, 0 } },
    { "eint",          FORT_DERIENG, { 1, 1, 1, 1, 1, 0 } },
    { "kint",          FORT_DERKENG, { 1, 1, 1, 1, 0, 0 } },
    { "vorticity",     FORT_DERVORT, { 1, 1, 1, 1, 0, 0 } },
    { "concentration", FORT_DERCON,  { 1, 0, 0, 0, 0, 1 } },
    { "dilitation",    FORT_DERDILI, { 1, 1, 1, 1, 0, 0 } },
    { "density",       FORT_DERDEN,  { 1, 0, 0, 0, 0, 0 } },
    { "mach",          FORT_DERMACH, { 1, 1, 1, 1, 1, 0 } },
    { "tracer",        FORT_DERTRAC, { 0, 0, 0, 0, 0, 1 } }
};



#endif /* BL_DERIVE_HCHEM */

#endif /*BL_DERIVE_IAMR*/

static const int NumDeriveFunc = sizeof(DerivedQuants)/sizeof(DeriveRecord);

#endif /*BL_SPACEDIM==3*/

#endif
