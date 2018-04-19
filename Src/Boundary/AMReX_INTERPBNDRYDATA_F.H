#ifndef AMREX_INTERPBNDRY_F_H_
#define AMREX_INTERPBNDRY_F_H_

#include <AMReX_BLFort.H>

#ifndef AMREX_XSDK
namespace amrex {
#endif
typedef void BDInterpFunc(amrex_real* bndr, AMREX_ARLIM_P(blo), AMREX_ARLIM_P(bhi),
                          const int* lo, const int* hi,
                          AMREX_ARLIM_P(cblo), AMREX_ARLIM_P(cbhi),
                          const int* nvar, const int* ratio,
                          const int* not_covered,
                          const int* mask,  AMREX_ARLIM_P(mlo), AMREX_ARLIM_P(mhi),
                          const amrex_real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
                          amrex_real* derives, int* max_order);
#ifndef AMREX_XSDK
}
#endif

#ifdef __cplusplus
extern "C"
{
#endif
    void amrex_bdinterpxlo(amrex_real* bndr, AMREX_ARLIM_P(blo), AMREX_ARLIM_P(bhi),
                          const int* lo, const int* hi,
                          AMREX_ARLIM_P(cblo), AMREX_ARLIM_P(cbhi),
                          const int* nvar, const int* ratio,
                          const int* not_covered,
                          const int* mask,  AMREX_ARLIM_P(mlo), AMREX_ARLIM_P(mhi),
                          const amrex_real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
                          amrex_real* derives, int* max_order);
    void amrex_bdinterpxhi(amrex_real* bndr, AMREX_ARLIM_P(blo), AMREX_ARLIM_P(bhi),
                          const int* lo, const int* hi,
                          AMREX_ARLIM_P(cblo), AMREX_ARLIM_P(cbhi),
                          const int* nvar, const int* ratio,
                          const int* not_covered,
                          const int* mask,  AMREX_ARLIM_P(mlo), AMREX_ARLIM_P(mhi),
                          const amrex_real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
                          amrex_real* derives, int* max_order);

#if (AMREX_SPACEDIM > 1)    
    void amrex_bdinterpylo(amrex_real* bndr, AMREX_ARLIM_P(blo), AMREX_ARLIM_P(bhi),
                          const int* lo, const int* hi,
                          AMREX_ARLIM_P(cblo), AMREX_ARLIM_P(cbhi),
                          const int* nvar, const int* ratio,
                          const int* not_covered,
                          const int* mask,  AMREX_ARLIM_P(mlo), AMREX_ARLIM_P(mhi),
                          const amrex_real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
                          amrex_real* derives, int* max_order);
    void amrex_bdinterpyhi(amrex_real* bndr, AMREX_ARLIM_P(blo), AMREX_ARLIM_P(bhi),
                          const int* lo, const int* hi,
                          AMREX_ARLIM_P(cblo), AMREX_ARLIM_P(cbhi),
                          const int* nvar, const int* ratio,
                          const int* not_covered,
                          const int* mask,  AMREX_ARLIM_P(mlo), AMREX_ARLIM_P(mhi),
                          const amrex_real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
                          amrex_real* derives, int* max_order);
#endif    

#if (AMREX_SPACEDIM > 2)    
    void amrex_bdinterpzlo(amrex_real* bndr, AMREX_ARLIM_P(blo), AMREX_ARLIM_P(bhi),
                          const int* lo, const int* hi,
                          AMREX_ARLIM_P(cblo), AMREX_ARLIM_P(cbhi),
                          const int* nvar, const int* ratio,
                          const int* not_covered,
                          const int* mask,  AMREX_ARLIM_P(mlo), AMREX_ARLIM_P(mhi),
                          const amrex_real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
                          amrex_real* derives, int* max_order);
    void amrex_bdinterpzhi(amrex_real* bndr, AMREX_ARLIM_P(blo), AMREX_ARLIM_P(bhi),
                          const int* lo, const int* hi,
                          AMREX_ARLIM_P(cblo), AMREX_ARLIM_P(cbhi),
                          const int* nvar, const int* ratio,
                          const int* not_covered,
                          const int* mask,  AMREX_ARLIM_P(mlo), AMREX_ARLIM_P(mhi),
                          const amrex_real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
                          amrex_real* derives, int* max_order);
#endif    

#ifdef __cplusplus
}
#endif

#endif /*_INTERPBNDRY_F_H_*/
