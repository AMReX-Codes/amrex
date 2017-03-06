#ifndef MYFUNC_F_H__ 
#define MUFUNC_F_H__ 

#include <AMReX_BLFort.H>

extern "C"
{
  void init_phi(amrex::Real* data, const int* lo, const int* hi, const int* ng, 
		const amrex::Real* dx, const amrex::Real* prob_lo, const amrex::Real* prob_hi);

  void compute_flux(amrex::Real* phi, const int* ng_p,
		    amrex::Real* fluxx, 
		    amrex::Real* fluxy,
#if (BL_SPACEDIM == 3)   
		    amrex::Real* fluxz,
#endif
		    const int* ng_f, const int* lo, const int* hi, 
		    const int* dlo, const int* dhi, 
		    const int* bc, const amrex::Real* dx);
  
  void update_phi(amrex::Real* phiold, amrex::Real* phinew, const int* ng_p,
		  amrex::Real* fluxx, 
		  amrex::Real* fluxy,
#if (BL_SPACEDIM == 3)   
		  amrex::Real* fluxz,
#endif
		  const int* ng_f, const int* lo, const int* hi, const amrex::Real* dx, const amrex::Real* dt);

    void phifill (BL_FORT_FAB_ARG_3D(state),
		  const int* dlo, const int* dhi,
		  const amrex::Real* dx, const amrex::Real* glo, 
		  const amrex::Real* time, const int* bc);
}
#endif
