#ifndef MYFUNC_F_H__ 
#define MUFUNC_F_H__ 

#include <BLFort.H>

extern "C"
{
  void init_phi(Real* data, const int* lo, const int* hi, const int* ng, 
		const Real* dx, const Real* prob_lo, const Real* prob_hi);

  void compute_flux(Real* phi, const int* ng_p,
		    Real* fluxx, 
		    Real* fluxy,
#if (BL_SPACEDIM == 3)   
		    Real* fluxz,
#endif
		    const int* ng_f, const int* lo, const int* hi, 
		    const int* dlo, const int* dhi, 
		    const int* bc, const Real* dx);
  
  void update_phi(Real* phiold, Real* phinew, const int* ng_p,
		  Real* fluxx, 
		  Real* fluxy,
#if (BL_SPACEDIM == 3)   
		  Real* fluxz,
#endif
		  const int* ng_f, const int* lo, const int* hi, const Real* dx, const Real* dt);

    void phifill (BL_FORT_FAB_ARG_3D(state),
		  const int* dlo, const int* dhi,
		  const Real* dx, const Real* glo, 
		  const Real* time, const int* bc);
}
#endif
