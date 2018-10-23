#ifndef AMREX_PROB_AMR_F_H_
#define AMREX_PROB_AMR_F_H_

#include <AMReX_REAL.H>

#ifdef __cplusplus
extern "C"
{
#endif
     void amrex_probinit (const int* init,
			  const int* name,
			  const int* namelen,
			  const amrex_real* problo,
			  const amrex_real* probhi);
#ifdef __cplusplus
}
#endif

#endif /*_PROB_AMR_F_H_*/
