#ifndef AMREX_INTERP_NAMES_H
#define AMREX_INTERP_NAMES_H

#ifdef BL_LANG_FORT
#    define FORT_NBINTERP    nbinterp
#    define FORT_CBINTERP    cbinterp
#    define FORT_CCINTERP    ccinterp
#    define FORT_LINCCINTERP linccinterp
#    define FORT_CQINTERP    cqinterp
#    define FORT_CCINTERP2   ccinterp2
#    define FORT_PCINTERP    pcinterp
#    define FORT_PROTECT_INTERP    printerp
#    define FORT_QUARTINTERP quartinterp
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_NBINTERP    NBINTERP
#    define FORT_CBINTERP    CBINTERP
#    define FORT_CCINTERP    CCINTERP
#    define FORT_LINCCINTERP LINCCINTERP
#    define FORT_CQINTERP    CQINTERP
#    define FORT_CCINTERP2   CCINTERP2
#    define FORT_PCINTERP    PCINTERP
#    define FORT_PROTECT_INTERP    PRINTERP
#    define FORT_QUARTINTERP QUARTINTERP
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_NBINTERP    nbinterp
#    define FORT_CBINTERP    cbinterp
#    define FORT_CCINTERP    ccinterp
#    define FORT_LINCCINTERP linccinterp
#    define FORT_CQINTERP    cqinterp
#    define FORT_CCINTERP2   ccinterp2
#    define FORT_PCINTERP    pcinterp
#    define FORT_PROTECT_INTERP    printerp
#    define FORT_QUARTINTERP quartinterp
#  else
#    define FORT_NBINTERP    nbinterp_
#    define FORT_CBINTERP    cbinterp_
#    define FORT_CCINTERP    ccinterp_
#    define FORT_LINCCINTERP linccinterp_
#    define FORT_CQINTERP    cqinterp_
#    define FORT_CCINTERP2   ccinterp2_
#    define FORT_PCINTERP    pcinterp_
#    define FORT_PROTECT_INTERP    printerp_
#    define FORT_QUARTINTERP quartinterp_
#  endif
#  endif

#include <AMReX_ArrayLim.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void FORT_NBINTERP (const amrex_real* crse, ARLIM_P(clo), ARLIM_P(chi),
                        ARLIM_P(cblo), ARLIM_P(cbhi),
                        amrex_real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                        ARLIM_P(fblo), ARLIM_P(fbhi),
                        AMREX_D_DECL(const int* lrx,const int* lry,const int* lrz),
                        const int* nvar, amrex_real* slope, const int* num_slope,
                        const int* actual_comp, const int* actual_state);
      
    void FORT_CBINTERP (const amrex_real* crse, ARLIM_P(clo), ARLIM_P(chi),
                        ARLIM_P(cblo), ARLIM_P(cbhi),
                        amrex_real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                        ARLIM_P(fblo), ARLIM_P(fbhi),
                        AMREX_D_DECL(const int* lrx,const int* lry,const int* lrz),
                        const int* nvar, amrex_real* slope, const int* num_slope,
                        amrex_real* strip, const int* strip_lo, const int* strip_hi,
                        const int* actual_comp, const int* actual_state);

    void FORT_CCINTERP (amrex_real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                        ARLIM_P(fblo), ARLIM_P(fbhi),
                        const int* nvar,
                        AMREX_D_DECL(const int* lrx,const int* lry,const int* lrz),
                        const amrex_real* crse, const int* clo, const int* chi,
                        ARLIM_P(cblo), ARLIM_P(cbhi),
                        const int* fslo, const int* fshi,
                        amrex_real* cslope, const int* clen,
                        amrex_real* fslope, amrex_real* fdat, const int* flen,
                        amrex_real* voff,
                        const int* bc, const int* lim_slope,
                        AMREX_D_DECL(const amrex_real* fvcx,const amrex_real* fvcy, const amrex_real* fvcz),
                        AMREX_D_DECL(const amrex_real* cvcx,const amrex_real* cvcy, const amrex_real* cvcz),
                        amrex_real* cmax, amrex_real* cmin, amrex_real* alpha,
                        const int* actual_comp, const int* actual_state);

    void FORT_CCINTERP2 (amrex_real* fine, const int* flo, const int* fhi,
                         ARLIM_P(fblo), ARLIM_P(fbhi),
                         const int* nvar,
                         AMREX_D_DECL(const int* lrx,const int* lry,const int* lrz),
                         const amrex_real* crse, const int* clo, const int* chi,
                         ARLIM_P(cblo), ARLIM_P(cbhi),
                         const int* fslo, const int* fshi,
                         amrex_real* cslope, const int* clen,
                         amrex_real* fslope, amrex_real* fdat, const int* flen,
                         amrex_real* voff,
                         const int* bc, const int* lim_slope,
                         AMREX_D_DECL(const amrex_real* fvcx,const amrex_real* fvcy, const amrex_real* fvcz),
                         AMREX_D_DECL(const amrex_real* cvcx,const amrex_real* cvcy, const amrex_real* cvcz),
                         const int* actual_comp, const int* actual_state);

    void FORT_LINCCINTERP (amrex_real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                           const int* fblo, const int* fbhi,
                           ARLIM_P(fvcblo), ARLIM_P(fvcbhi),
                           const amrex_real* crse, ARLIM_P(clo), ARLIM_P(chi),
                           ARLIM_P(cvcblo), ARLIM_P(cvcbhi),
                           amrex_real* ucc_xslopes, amrex_real* lcc_xslopes,
                           amrex_real* xslope_factors,
#if (AMREX_SPACEDIM>=2)
                           amrex_real* ucc_yslopes, amrex_real* lcc_yslopes,
                           amrex_real* yslope_factors,
#endif
#if (AMREX_SPACEDIM==3)
                           amrex_real* ucc_zslopes, amrex_real* lcc_zslopes,
                           amrex_real* zslope_factors,
#endif
                           ARLIM_P(csblo), ARLIM_P(csbhi),
                           const int* csblo, const int* csbhi,
                           const int* nvar,
                           AMREX_D_DECL(const int* lrx,const int* lry,const int* lrz),
                           const int* bc, 
                           const int* lim_slope, const int* lim_limit,
                           AMREX_D_DECL(const amrex_real* fvcx,const amrex_real* fvcy, const amrex_real* fvcz),
                           AMREX_D_DECL(const amrex_real* cvcx,const amrex_real* cvcy, const amrex_real* cvcz),
                           AMREX_D_DECL(const amrex_real* voffx,const amrex_real* voffy, const amrex_real* voffz),
                           amrex_real* alpha, amrex_real* cmax, amrex_real* cmin,
                           const int* actual_comp, const int* actual_state);

    void FORT_CQINTERP (amrex_real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                        ARLIM_P(fblo), ARLIM_P(fbhi),
                        const int* nvar,
                        AMREX_D_DECL(const int* lrx,const int* lry,const int* lrz),
                        const amrex_real* crse, const int* clo, const int* chi,
                        ARLIM_P(cblo), ARLIM_P(cbhi),
                        const int* fslo, const int* fshi,
                        amrex_real* cslope, const int* clen,
                        amrex_real* fslope, amrex_real* fdat, const int* flen,
                        amrex_real* voff,
                        const int* bc, const int* lim_slope,
                        AMREX_D_DECL(const amrex_real* fvcx,const amrex_real* fvcy, const amrex_real* fvcz),
                        AMREX_D_DECL(const amrex_real* cvcx,const amrex_real* cvcy, const amrex_real* cvcz),
                        const int* actual_comp, const int* actual_state);

    void FORT_PCINTERP (const amrex_real* crse, ARLIM_P(clo), ARLIM_P(chi),
                        const int* cblo, const int* cbhi,
                        amrex_real*  fine, ARLIM_P(flo), ARLIM_P(fhi),
                        const int* fblo, const int* fbhi,
                        const int* longdir,
                        AMREX_D_DECL(const int* lrx,const int* lry,const int* lrz),
                        const int* nvar, amrex_real* strip,
                        const int* strip_lo, const int* strip_hi,
                        const int* actual_comp, const int* actual_state);

    void FORT_PROTECT_INTERP (      amrex_real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                              const  int* fblo, const int* fbhi,
                              const amrex_real* crse, ARLIM_P(clo), ARLIM_P(chi),
                              const  int* cblo, const int* cbhi,
#if (AMREX_SPACEDIM == 2)
                              const amrex_real* fvcx, const amrex_real* fvcy,
                              ARLIM_P(fblo), ARLIM_P(fbhi),
                              const amrex_real* cvcx, const amrex_real* cvcy,
                              ARLIM_P(cblo), ARLIM_P(cbhi),
#endif
                              amrex_real* state, ARLIM_P(slo), ARLIM_P(shi),
                              const int* nvar,
                              AMREX_D_DECL(const int* lrx,const int* lry,const int* lrz),
                              const int* bc); 

    void FORT_QUARTINTERP (amrex_real* fine, ARLIM_P(flo), ARLIM_P(fhi),
			   const int* fblo, const int* fbhi, 
			   const int* fb2lo, const int* fb2hi,
			   const amrex_real* crse, ARLIM_P(clo), ARLIM_P(chi),
			   const int* cblo, const int* cbhi, 
			   const int* cb2lo, const int* cb2hi,
			   const int* nvar, 
			   AMREX_D_DECL(const int* lrx,const int* lry,const int* lrz),
			   AMREX_D_DECL(amrex_real* ftmp, amrex_real* ctmp, amrex_real* ctmp2),
			   const int* bc, const int* actual_comp, const int* actual_state);
#ifdef __cplusplus
  }
#endif
#endif

#endif /*_INTERP_NAMES_H*/
