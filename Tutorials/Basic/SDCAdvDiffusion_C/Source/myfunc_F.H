#ifndef MYFUNC_F_H__ 
#define MYFUNC_F_H__ 

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void init_phi(const int* lo, const int* hi,
		amrex_real* data, const int* dlo, const int* dhi,
		const amrex_real* dx, const amrex_real* prob_lo, const amrex_real* prob_hi);
  void pf_quadrature(int* qtype_in, int* nnodes, int* nnodes0,
		     amrex_real* nodes, int* nflags,amrex_real* qmats);
  
  void err_phi(const int* lo, const int* hi,
	       amrex_real* data, const int* dlo, const int* dhi,
	       const amrex_real* dx, const amrex_real* prob_lo, const amrex_real* prob_hi,
	       amrex_real* time,amrex_real* v,amrex_real* nu);
  
  void evaluate_f(const int* lo, const int* hi,
		 const int* domlo, const int* domhi,
		 amrex_real* phi, const int* philo, const int* phihi,
		 amrex_real* fluxx, const int* fxlo, const int* fxhi,
		 amrex_real* fluxy, const int* fylo, const int* fyhi,
#if (BL_SPACEDIM == 3)   
		 amrex_real* fluxz, const int* fzlo, const int* fzhi,
#endif
		 amrex_real* f, const int* flo, const int* fhi,
		 const amrex_real* dx,amrex_real* v,amrex_real* nu,int *npiece);
  
  
#ifdef __cplusplus
}
#endif

#endif
