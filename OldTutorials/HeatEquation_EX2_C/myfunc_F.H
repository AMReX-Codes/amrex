#ifndef MYFUNC_F_H__ 
#define MUFUNC_F_H__ 

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void init_phi(amrex_real* data, const int* lo, const int* hi, const int* ng, 
		const amrex_real* dx, const amrex_real* prob_lo, const amrex_real* prob_hi);

  void compute_flux(amrex_real* phi, const int* ng_p,
		    amrex_real* fluxx, 
		    amrex_real* fluxy,
#if (BL_SPACEDIM == 3)   
		    amrex_real* fluxz,
#endif
		    const int* ng_f, const int* lo, const int* hi, 
		    const int* dlo, const int* dhi, 
		    const int* bc, const amrex_real* dx);
  
  void update_phi(amrex_real* phiold, amrex_real* phinew, const int* ng_p,
		  amrex_real* fluxx, 
		  amrex_real* fluxy,
#if (BL_SPACEDIM == 3)   
		  amrex_real* fluxz,
#endif
		  const int* ng_f, const int* lo, const int* hi, const amrex_real* dx, const amrex_real* dt);

    void phifill (BL_FORT_FAB_ARG_3D(state),
		  const int* dlo, const int* dhi,
		  const amrex_real* dx, const amrex_real* glo, 
		  const amrex_real* time, const int* bc);

#ifdef __cplusplus
}
#endif

#endif
