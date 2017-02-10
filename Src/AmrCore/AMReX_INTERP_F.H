#ifndef _INTERP_NAMES_H
#define _INTERP_NAMES_H

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

extern "C"
{
    void FORT_NBINTERP (const amrex::Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                        ARLIM_P(cblo), ARLIM_P(cbhi),
                        amrex::Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                        ARLIM_P(fblo), ARLIM_P(fbhi),
                        D_DECL(const int* lrx,const int* lry,const int* lrz),
                        const int* nvar, amrex::Real* slope, const int* num_slope,
                        const int* actual_comp, const int* actual_state);
      
    void FORT_CBINTERP (const amrex::Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                        ARLIM_P(cblo), ARLIM_P(cbhi),
                        amrex::Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                        ARLIM_P(fblo), ARLIM_P(fbhi),
                        D_DECL(const int* lrx,const int* lry,const int* lrz),
                        const int* nvar, amrex::Real* slope, const int* num_slope,
                        amrex::Real* strip, const int* strip_lo, const int* strip_hi,
                        const int* actual_comp, const int* actual_state);

    void FORT_CCINTERP (amrex::Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                        ARLIM_P(fblo), ARLIM_P(fbhi),
                        const int* nvar,
                        D_DECL(const int* lrx,const int* lry,const int* lrz),
                        const amrex::Real* crse, const int* clo, const int* chi,
                        ARLIM_P(cblo), ARLIM_P(cbhi),
                        const int* fslo, const int* fshi,
                        amrex::Real* cslope, const int* clen,
                        amrex::Real* fslope, amrex::Real* fdat, const int* flen,
                        amrex::Real* voff,
                        const int* bc, const int* lim_slope,
                        D_DECL(const amrex::Real* fvcx,const amrex::Real* fvcy, const amrex::Real* fvcz),
                        D_DECL(const amrex::Real* cvcx,const amrex::Real* cvcy, const amrex::Real* cvcz),
                        amrex::Real* cmax, amrex::Real* cmin, amrex::Real* alpha,
                        const int* actual_comp, const int* actual_state);

    void FORT_CCINTERP2 (amrex::Real* fine, const int* flo, const int* fhi,
                         ARLIM_P(fblo), ARLIM_P(fbhi),
                         const int* nvar,
                         D_DECL(const int* lrx,const int* lry,const int* lrz),
                         const amrex::Real* crse, const int* clo, const int* chi,
                         ARLIM_P(cblo), ARLIM_P(cbhi),
                         const int* fslo, const int* fshi,
                         amrex::Real* cslope, const int* clen,
                         amrex::Real* fslope, amrex::Real* fdat, const int* flen,
                         amrex::Real* voff,
                         const int* bc, const int* lim_slope,
                         D_DECL(const amrex::Real* fvcx,const amrex::Real* fvcy, const amrex::Real* fvcz),
                         D_DECL(const amrex::Real* cvcx,const amrex::Real* cvcy, const amrex::Real* cvcz),
                         const int* actual_comp, const int* actual_state);

    void FORT_LINCCINTERP (amrex::Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                           const int* fblo, const int* fbhi,
                           ARLIM_P(fvcblo), ARLIM_P(fvcbhi),
                           const amrex::Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                           ARLIM_P(cvcblo), ARLIM_P(cvcbhi),
                           amrex::Real* ucc_xslopes, amrex::Real* lcc_xslopes,
                           amrex::Real* xslope_factors,
#if (BL_SPACEDIM>=2)
                           amrex::Real* ucc_yslopes, amrex::Real* lcc_yslopes,
                           amrex::Real* yslope_factors,
#endif
#if (BL_SPACEDIM==3)
                           amrex::Real* ucc_zslopes, amrex::Real* lcc_zslopes,
                           amrex::Real* zslope_factors,
#endif
                           ARLIM_P(csblo), ARLIM_P(csbhi),
                           const int* csblo, const int* csbhi,
                           const int* nvar,
                           D_DECL(const int* lrx,const int* lry,const int* lrz),
                           const int* bc, 
                           const int* lim_slope, const int* lim_limit,
                           D_DECL(const amrex::Real* fvcx,const amrex::Real* fvcy, const amrex::Real* fvcz),
                           D_DECL(const amrex::Real* cvcx,const amrex::Real* cvcy, const amrex::Real* cvcz),
                           D_DECL(const amrex::Real* voffx,const amrex::Real* voffy, const amrex::Real* voffz),
                           amrex::Real* alpha, amrex::Real* cmax, amrex::Real* cmin,
                           const int* actual_comp, const int* actual_state);

    void FORT_CQINTERP (amrex::Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                        ARLIM_P(fblo), ARLIM_P(fbhi),
                        const int* nvar,
                        D_DECL(const int* lrx,const int* lry,const int* lrz),
                        const amrex::Real* crse, const int* clo, const int* chi,
                        ARLIM_P(cblo), ARLIM_P(cbhi),
                        const int* fslo, const int* fshi,
                        amrex::Real* cslope, const int* clen,
                        amrex::Real* fslope, amrex::Real* fdat, const int* flen,
                        amrex::Real* voff,
                        const int* bc, const int* lim_slope,
                        D_DECL(const amrex::Real* fvcx,const amrex::Real* fvcy, const amrex::Real* fvcz),
                        D_DECL(const amrex::Real* cvcx,const amrex::Real* cvcy, const amrex::Real* cvcz),
                        const int* actual_comp, const int* actual_state);

    void FORT_PCINTERP (const amrex::Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                        const int* cblo, const int* cbhi,
                        amrex::Real*  fine, ARLIM_P(flo), ARLIM_P(fhi),
                        const int* fblo, const int* fbhi,
                        const int* longdir,
                        D_DECL(const int* lrx,const int* lry,const int* lrz),
                        const int* nvar, amrex::Real* strip,
                        const int* strip_lo, const int* strip_hi,
                        const int* actual_comp, const int* actual_state);

    void FORT_PROTECT_INTERP (      amrex::Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                              const  int* fblo, const int* fbhi,
                              const amrex::Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                              const  int* cblo, const int* cbhi,
#if (BL_SPACEDIM == 2)
                              const amrex::Real* fvcx, const amrex::Real* fvcy,
                              ARLIM_P(fblo), ARLIM_P(fbhi),
                              const amrex::Real* cvcx, const amrex::Real* cvcy,
                              ARLIM_P(cblo), ARLIM_P(cbhi),
#endif
                              amrex::Real* state, ARLIM_P(slo), ARLIM_P(shi),
                              const int* nvar,
                              D_DECL(const int* lrx,const int* lry,const int* lrz),
                              const int* bc); 

    void FORT_QUARTINTERP (amrex::Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
			   const int* fblo, const int* fbhi, 
			   const int* fb2lo, const int* fb2hi,
			   const amrex::Real* crse, ARLIM_P(clo), ARLIM_P(chi),
			   const int* cblo, const int* cbhi, 
			   const int* cb2lo, const int* cb2hi,
			   const int* nvar, 
			   D_DECL(const int* lrx,const int* lry,const int* lrz),
			   D_DECL(amrex::Real* ftmp, amrex::Real* ctmp, amrex::Real* ctmp2),
			   const int* bc, const int* actual_comp, const int* actual_state);
  }
#endif

#endif /*_INTERP_NAMES_H*/
