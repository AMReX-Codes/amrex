#ifndef AMREX_INTERPBNDRY_F_H_
#define AMREX_INTERPBNDRY_F_H_

#include <AMReX_BLFort.H>

typedef void BDInterpFunc(amrex_real* bndr, ARLIM_P(blo), ARLIM_P(bhi),
                          const int* lo, const int* hi,
                          ARLIM_P(cblo), ARLIM_P(cbhi),
                          const int* nvar, const int* ratio,
                          const int* not_covered,
                          const int* mask,  ARLIM_P(mlo), ARLIM_P(mhi),
                          const amrex_real* crse, ARLIM_P(clo), ARLIM_P(chi),
                          amrex_real* derives, int* max_order);

#ifdef __cplusplus
extern "C"
{
#endif
    BDInterpFunc amrex_bdinterpxlo;
    BDInterpFunc amrex_bdinterpxhi;

#if (AMREX_SPACEDIM > 1)    
    BDInterpFunc amrex_bdinterpylo;
    BDInterpFunc amrex_bdinterpyhi;
#endif    

#if (AMREX_SPACEDIM > 2)    
    BDInterpFunc amrex_bdinterpzlo;
    BDInterpFunc amrex_bdinterpzhi;
#endif    

#ifdef __cplusplus
}
#endif

#endif /*_INTERPBNDRY_F_H_*/
