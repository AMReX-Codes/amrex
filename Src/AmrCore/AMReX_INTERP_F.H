#ifndef AMREX_INTERP_NAMES_H
#define AMREX_INTERP_NAMES_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void amrex_cbinterp (const amrex_real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
                        AMREX_ARLIM_P(cblo), AMREX_ARLIM_P(cbhi),
                        amrex_real* fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
                        AMREX_ARLIM_P(fblo), AMREX_ARLIM_P(fbhi),
                        AMREX_D_DECL(const int* lrx,const int* lry,const int* lrz),
                        const int* nvar, amrex_real* slope, const int* num_slope,
                        amrex_real* strip, const int* strip_lo, const int* strip_hi,
                        const int* actual_comp, const int* actual_state);

    void amrex_cqinterp (amrex_real* fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
                        AMREX_ARLIM_P(fblo), AMREX_ARLIM_P(fbhi),
                        const int* nvar,
                        AMREX_D_DECL(const int* lrx,const int* lry,const int* lrz),
                        const amrex_real* crse, const int* clo, const int* chi,
                        AMREX_ARLIM_P(cblo), AMREX_ARLIM_P(cbhi),
                        const int* fslo, const int* fshi,
                        amrex_real* cslope, const int* clen,
                        amrex_real* fslope, amrex_real* fdat, const int* flen,
                        amrex_real* voff,
                        const int* bc, const int* lim_slope,
                        AMREX_D_DECL(const amrex_real* fvcx,const amrex_real* fvcy, const amrex_real* fvcz),
                        AMREX_D_DECL(const amrex_real* cvcx,const amrex_real* cvcy, const amrex_real* cvcz),
                        const int* actual_comp, const int* actual_state);

    void amrex_protect_interp (      amrex_real* fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
                              const  int* fblo, const int* fbhi,
                              const amrex_real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
                              const  int* cblo, const int* cbhi,
#if (AMREX_SPACEDIM == 2)
                              const amrex_real* fvcx, const amrex_real* fvcy,
                              AMREX_ARLIM_P(fblo), AMREX_ARLIM_P(fbhi),
                              const amrex_real* cvcx, const amrex_real* cvcy,
                              AMREX_ARLIM_P(cblo), AMREX_ARLIM_P(cbhi),
#endif
                              amrex_real* state, AMREX_ARLIM_P(slo), AMREX_ARLIM_P(shi),
                              const int* nvar,
                              AMREX_D_DECL(const int* lrx,const int* lry,const int* lrz),
                              const int* bc); 

    void amrex_quartinterp (amrex_real* fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
			   const int* fblo, const int* fbhi, 
			   const int* fb2lo, const int* fb2hi,
			   const amrex_real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
			   const int* cblo, const int* cbhi, 
			   const int* cb2lo, const int* cb2hi,
			   const int* nvar, 
			   AMREX_D_DECL(const int* lrx,const int* lry,const int* lrz),
			   AMREX_D_DECL(amrex_real* ftmp, amrex_real* ctmp, amrex_real* ctmp2),
			   const int* bc, const int* actual_comp, const int* actual_state);
#ifdef __cplusplus
  }
#endif


#endif /*_INTERP_NAMES_H*/
