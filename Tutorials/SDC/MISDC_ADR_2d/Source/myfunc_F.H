#ifndef MYFUNC_F_H_
#define MYFUNC_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void init_phi(const int* lo, const int* hi,
		amrex_real* data, const int* dlo, const int* dhi,
		const amrex_real* dx, const amrex_real* prob_lo,
		const amrex_real* prob_hi);
  
  void err_phi(const int* lo, const int* hi,
	       amrex_real* data, const int* dlo, const int* dhi,
	       const amrex_real* dx, const amrex_real* prob_lo, const amrex_real* prob_hi,
	       amrex_real* a,amrex_real* d,amrex_real* r,amrex_real* time);

  void SDC_feval_F(const int* lo, const int* hi,
		   const int* domlo, const int* domhi,
		   amrex_real* phi, const int* philo, const int* phihi,
		   amrex_real* fluxx, const int* fxlo, const int* fxhi,
		   amrex_real* fluxy, const int* fylo, const int* fyhi,
#if (BL_SPACEDIM == 3)   
		   amrex_real* fluxz, const int* fzlo, const int* fzhi,
#endif
		   amrex_real* f, const int* flo, const int* fhi,
		   const amrex_real* dx,const amrex_real* a,
		   const amrex_real* d,const amrex_real* r, int *n);
  void SDC_fcomp_reaction_F(const int* lo, const int* hi,
                      const int* domlo, const int* domhi,
                      amrex_real* phi, const int* philo, const int* phihi,
                      amrex_real* rhs, const int* rhslo, const int* rhshi,			    
                      amrex_real* f, const int* flo, const int* fhi,
   		    const amrex_real* dtq, int *n);
    
#ifdef __cplusplus
}
#endif

#endif
