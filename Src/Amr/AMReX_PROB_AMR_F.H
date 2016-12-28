#ifndef _PROB_AMR_F_H_
#define _PROB_AMR_F_H_

extern "C"
{
     void amrex_probinit (const int* init,
			  const int* name,
			  const int* namelen,
			  const Real* problo,
			  const Real* probhi);
}

#endif /*_PROB_AMR_F_H_*/
