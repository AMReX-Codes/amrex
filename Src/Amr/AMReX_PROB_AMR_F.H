#ifndef _PROB_AMR_F_H_
#define _PROB_AMR_F_H_

#include <AMReX_REAL.H>

#ifdef __cplusplus
extern "C"
{
#endif
     void amrex_probinit (const int* init,
			  const int* name,
			  const int* namelen,
			  const amrex::Real* problo,
			  const amrex::Real* probhi);
#ifdef __cplusplus
}
#endif

#endif /*_PROB_AMR_F_H_*/
