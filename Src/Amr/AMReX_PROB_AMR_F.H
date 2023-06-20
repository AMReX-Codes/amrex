#ifndef AMREX_PROB_AMR_F_H_
#define AMREX_PROB_AMR_F_H_
#include <AMReX_Config.H>

#include <AMReX_REAL.H>
#include <AMReX_Extension.H>

#ifdef __cplusplus
extern "C"
{
#endif
     void amrex_probinit (const int* init,
                          const int* name,
                          const int* namelen,
                          const amrex_real* problo,
                          const amrex_real* probhi) AMREX_ATTRIBUTE_WEAK;
#ifdef __cplusplus
}
#endif

#endif /*_PROB_AMR_F_H_*/
